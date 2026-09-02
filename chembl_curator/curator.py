# chembl_curator/curator.py

import sqlite3
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass

from .downloader import ChEMBLDownloader
from .filters import ActivityFilter, CompoundFilter
from .config import CurationConfig
from .labeler import Label, Measurement, classify_all
from .utils import setup_logging


@dataclass
class CurationResults:
    total_activities: int
    filtered_activities: int
    total_proteins: int
    total_compounds: int
    output_directory: Path
    total_inactives: int = 0
    total_conflicts: int = 0
    total_measured: int = 0


class ChEMBLCurator:    
    def __init__(
        self, 
        config: Optional[CurationConfig] = None,
        log_level: str = "INFO"
    ):
        """
        Args:
            config: Curation configs like activity thresholds
            log_level: DEBUG, INFO, WARNING, ERROR for logger
        """
        self.config = config or CurationConfig()
        self.logger = setup_logging(log_level)
        self.downloader = ChEMBLDownloader()
        self.activity_filter = ActivityFilter(self.config)
        self.compound_filter = CompoundFilter(self.config)
        
    def download_database(
        self, 
        output_dir: Optional[Path] = None,
        force_download: bool = False
    ) -> Path:
        """
        Args:
            output_dir: Directory to save database
            force_download: Force re-download even if exists
            
        Returns:
            Path to database file
        """
        return self.downloader.download_sqlite(
            output_dir=output_dir,
            force_download=force_download
        )
        
    def run_pipeline(
        self, 
        database_path: Optional[Path] = None,
        output_dir: Path = Path("./curated_chembl")
    ) -> CurationResults:
        """
        Args:
            database_path: Path to ChEMBL SQLite database
            output_dir: Output directory for curated data
            
        Returns:
            CurationResults with statistics
        """
        self.logger.info("Start curation for ChEMBL database")
        
        # Download database if not provided
        if database_path is None:
            database_path = self.download_database()
            
        conn = sqlite3.connect(database_path)
        
        try:
            activities = self._extract_activities(conn)
            self.logger.info(f"Extracted {len(activities)} activity records")
            
            chembl_to_uniprot = self._get_uniprot_mapping(conn)
            chembl_to_smiles = self._get_smiles_mapping(conn)
            
            results, active_pairs = self._process_activities(
                activities,
                chembl_to_uniprot,
                chembl_to_smiles,
                output_dir
            )

            if self.config.extract_negatives:
                negatives = self._extract_negatives(conn)
                self.logger.info(f"Extracted {len(negatives)} negative-candidate records")
                self._process_negatives(
                    negatives,
                    chembl_to_uniprot,
                    chembl_to_smiles,
                    output_dir,
                    active_pairs,
                    results,
                )

            self.logger.info(f"Curation completed: {results.total_compounds} actives for {results.total_proteins} targets")
            return results
            
        finally:
            conn.close()
    
    def _shared_quality_conditions(self) -> List[str]:
        """Assay and record quality conditions applied to actives and negatives alike.

        Built in one place so the two extraction queries cannot drift apart.
        Deliberately excludes the pChEMBL floor (an active-only potency rule)
        and standard_type (see _extract_negatives).
        """
        conditions = []

        if self.config.require_standard_flag:
            conditions.append("a.standard_flag = 1")

        if self.config.exclude_invalid_data:
            conditions.append("a.data_validity_comment IS NULL")

        if self.config.exclude_duplicates:
            conditions.append("(a.potential_duplicate IS NULL OR a.potential_duplicate = 0)")

        if self.config.min_confidence_score is not None:
            conditions.append(f"ass.confidence_score >= {self.config.min_confidence_score}")

        if self.config.assay_types is not None:
            assay_type_list = ','.join(f"'{t}'" for t in self.config.assay_types)
            conditions.append(f"ass.assay_type IN ({assay_type_list})")

        if self.config.bao_formats is not None:
            bao_format_list = ','.join(f"'{b}'" for b in self.config.bao_formats)
            conditions.append(f"ass.bao_format IN ({bao_format_list})")

        return conditions

    def _extract_activities(self, conn: sqlite3.Connection) -> List[Tuple]:
        target_types = ','.join(f"'{t}'" for t in self.config.target_types)
        activity_types = ','.join(f"'{t}'" for t in self.config.activity_types)
        relations = ','.join(f"'{r}'" for r in self.config.relations)
        units = ','.join(f"'{u}'" for u in self.config.units)

        validity_conditions = self._shared_quality_conditions()

        # Add pChEMBL filter (active-only: negatives have no potency to threshold)
        if self.config.min_pchembl_value is not None:
            validity_conditions.append(f"a.pchembl_value >= {self.config.min_pchembl_value}")

        # Combine validity conditions
        validity_filter = ""
        if validity_conditions:
            validity_filter = "AND " + " AND ".join(validity_conditions)

        query = f"""
        SELECT
            td.chembl_id AS target_chembl_id,
            md.chembl_id AS compound_chembl_id,
            a.standard_type,
            a.standard_relation,
            a.standard_value,
            a.standard_units,
            a.pchembl_value
        FROM activities a
        JOIN molecule_dictionary md ON a.molregno = md.molregno
        JOIN assays ass ON a.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        WHERE
            td.target_type IN ({target_types})
            AND a.standard_type IN ({activity_types})
            AND a.standard_relation IN ({relations})
            AND a.standard_units IN ({units})
            AND a.standard_value IS NOT NULL
            {validity_filter}
        """
        return conn.execute(query).fetchall()
    
    def _extract_negatives(self, conn: sqlite3.Connection) -> List[Tuple]:
        """Pull every measurement usable for negative labelling.

        Uses the identical assay and record quality conditions as the active
        query. Three deliberate relaxations, each with a reason:

        - no pChEMBL floor and no relation filter, so censored (">") records
          and sub-threshold values arrive at all;
        - standard_value may be NULL, since an explicit inactive call often
          carries no number;
        - no standard_type restriction. An inactive compound has no IC50 to
          report, so its result is filed under "% Control", "Inhibition" or
          "Activity" instead. Restricting to potency types drops ~90% of the
          depositor inactive calls. labeler.concentration_nm() makes sure
          those non-potency numbers are never read as concentrations.
        """
        target_types = ','.join(f"'{t}'" for t in self.config.target_types)
        conditions = self._shared_quality_conditions()
        quality_filter = ("AND " + " AND ".join(conditions)) if conditions else ""

        query = f"""
        SELECT
            td.chembl_id AS target_chembl_id,
            md.chembl_id AS compound_chembl_id,
            a.standard_type,
            a.standard_relation,
            a.standard_value,
            a.standard_units,
            a.pchembl_value,
            a.activity_comment
        FROM activities a
        JOIN molecule_dictionary md ON a.molregno = md.molregno
        JOIN assays ass ON a.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        WHERE
            td.target_type IN ({target_types})
            {quality_filter}
        """
        return conn.execute(query).fetchall()

    def _get_uniprot_mapping(self, conn: sqlite3.Connection) -> Dict[str, str]:
        target_types = ','.join(f"'{t}'" for t in self.config.target_types)
        
        query = f"""
        SELECT DISTINCT 
            td.chembl_id AS target_chembl_id,
            cs.accession AS uniprot_id
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE 
            td.target_type IN ({target_types})
            AND cs.accession IS NOT NULL
        """
        results = conn.execute(query).fetchall()
        return {chembl_id: uniprot_id for chembl_id, uniprot_id in results}
    
    def _get_smiles_mapping(self, conn: sqlite3.Connection) -> Dict[str, str]:
        query = """
        SELECT 
            md.chembl_id,
            cs.canonical_smiles
        FROM molecule_dictionary md
        JOIN compound_structures cs ON md.molregno = cs.molregno
        WHERE cs.canonical_smiles IS NOT NULL
        """
        results = conn.execute(query).fetchall()
        return {chembl_id: smiles for chembl_id, smiles in results}
    
    def _process_activities(
        self,
        activities: List[Tuple],
        chembl_to_uniprot: Dict[str, str],
        chembl_to_smiles: Dict[str, str],
        output_dir: Path
    ) -> Tuple[CurationResults, Set[Tuple[str, str]]]:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # {uniprot: {chembl_id: (smiles, pchembl)}} — keeps best pChEMBL per compound per target
        protein_compounds = {}
        filtered_count = 0

        for activity in activities:
            target_chembl, compound_chembl, _std_type, _std_relation, std_value, std_units, pchembl = activity

            if not self.activity_filter.is_active(std_value, std_units):
                continue

            if target_chembl not in chembl_to_uniprot:
                continue
            uniprot_id = chembl_to_uniprot[target_chembl]

            if compound_chembl not in chembl_to_smiles:
                continue
            smiles = chembl_to_smiles[compound_chembl]

            if not self.compound_filter.is_valid(smiles):
                continue

            pchembl_f = float(pchembl) if pchembl is not None else 0.0

            if uniprot_id not in protein_compounds:
                protein_compounds[uniprot_id] = {}
            existing = protein_compounds[uniprot_id].get(compound_chembl)
            if existing is None or pchembl_f > existing[1]:
                protein_compounds[uniprot_id][compound_chembl] = (smiles, pchembl_f)
            filtered_count += 1

        total_compounds = self._write_output_files(protein_compounds, output_dir)

        active_pairs = {
            (uniprot_id, compound_id)
            for uniprot_id, compounds in protein_compounds.items()
            for compound_id in compounds
        }

        results = CurationResults(
            total_activities=len(activities),
            filtered_activities=filtered_count,
            total_proteins=len(protein_compounds),
            total_compounds=total_compounds,
            output_directory=output_dir
        )
        return results, active_pairs
    
    def _process_negatives(
        self,
        negatives: List[Tuple],
        chembl_to_uniprot: Dict[str, str],
        chembl_to_smiles: Dict[str, str],
        output_dir: Path,
        active_pairs: Set[Tuple[str, str]],
        results: CurationResults,
    ) -> None:
        """Label negative candidates and write inactives/measured/conflicts per target.

        A pair that the active pass accepted but that also carries inactive
        evidence is a contradiction: it is written to conflicts.tsv and kept
        out of both label sets.
        """
        measurements = [
            Measurement(
                target=chembl_to_uniprot[row[0]],
                compound=row[1],
                std_type=row[2],
                relation=row[3],
                value=row[4],
                units=row[5],
                pchembl=row[6],
                comment=row[7],
            )
            for row in negatives
            if row[0] in chembl_to_uniprot and row[1] in chembl_to_smiles
        ]
        self.logger.info(f"Labelling {len(measurements)} mapped negative records")

        verdicts = classify_all(
            measurements,
            active_max_nm=self.config.active_max_nm,
            inactive_min_nm=self.config.inactive_min_nm,
            conflict_decisive_nm=self.config.conflict_decisive_nm,
            potency_types=frozenset(self.config.activity_types),
        )

        inactives: Dict[str, List[Tuple[str, str]]] = {}
        conflicts: Dict[str, List[Tuple]] = {}
        measured: Dict[str, Set[str]] = {}
        n_contested = 0

        for (uniprot_id, compound_id), verdict in verdicts.items():
            measured.setdefault(uniprot_id, set()).add(compound_id)
            evidence = ";".join(verdict.evidence)
            if verdict.is_contested:
                n_contested += 1
            # A pair the active pass accepted must never also become an inactive.
            is_conflict = verdict.label is Label.CONFLICT or (
                verdict.is_inactive and (uniprot_id, compound_id) in active_pairs
            )
            if is_conflict:
                conflicts.setdefault(uniprot_id, []).append(
                    (compound_id, evidence, verdict.n_active_records,
                     verdict.n_inactive_records, verdict.best_active_nm)
                )
            elif verdict.is_inactive:
                inactives.setdefault(uniprot_id, []).append((compound_id, evidence))

        # Active pairs are measured too, even when this pass saw no negative row.
        for uniprot_id, compound_id in active_pairs:
            measured.setdefault(uniprot_id, set()).add(compound_id)

        conflict_ids = {
            (uniprot_id, row[0])
            for uniprot_id, rows in conflicts.items()
            for row in rows
        }

        for uniprot_id in measured:
            target_dir = output_dir / uniprot_id
            if not target_dir.exists():
                continue

            with open(target_dir / "inactives.tsv", 'w') as f:
                f.write("chembl_id\tevidence\tsmiles\n")
                for compound_id, evidence in sorted(inactives.get(uniprot_id, [])):
                    f.write(f"{compound_id}\t{evidence}\t{chembl_to_smiles[compound_id]}\n")

            with open(target_dir / "measured.tsv", 'w') as f:
                f.write("chembl_id\n")
                for compound_id in sorted(measured[uniprot_id]):
                    f.write(f"{compound_id}\n")

            rows = sorted(conflicts.get(uniprot_id, []))
            if rows:
                with open(target_dir / "conflicts.tsv", 'w') as f:
                    f.write("chembl_id\tevidence\tn_active\tn_inactive\tbest_active_nm\n")
                    for compound_id, evidence, n_act, n_inact, best in rows:
                        best_str = "" if best is None else f"{best:.4g}"
                        f.write(
                            f"{compound_id}\t{evidence}\t{n_act}\t{n_inact}\t{best_str}\n"
                        )

        results.total_inactives = sum(len(v) for v in inactives.values())
        results.total_conflicts = len(conflict_ids)
        results.total_measured = sum(len(v) for v in measured.values())
        self.logger.info(
            f"Negatives: {results.total_inactives} inactive pairs, "
            f"{results.total_conflicts} conflicts removed, "
            f"{n_contested} contested pairs kept active on potency, "
            f"{results.total_measured} measured pairs (decoy exclusion list)"
        )

    def _write_output_files(
        self,
        protein_compounds: Dict[str, Dict],
        output_dir: Path
    ) -> int:
        total_compounds = 0

        for uniprot_id, compounds in protein_compounds.items():
            target_dir = output_dir / uniprot_id
            smi_dir = target_dir / "comps" / "smiles"
            smi_dir.mkdir(parents=True, exist_ok=True)

            # Write individual .smi files (existing format)
            for compound_id, (smiles, _) in compounds.items():
                smi_file = smi_dir / f"{compound_id}.smi"
                with open(smi_file, 'w') as f:
                    f.write(smiles)

            # Write actives.tsv (chembl_id, pchembl, smiles) — used by downstream stages
            tsv_path = target_dir / "actives.tsv"
            with open(tsv_path, 'w') as f:
                f.write("chembl_id\tpchembl\tsmiles\n")
                for compound_id, (smiles, pchembl) in sorted(compounds.items()):
                    f.write(f"{compound_id}\t{pchembl:.4f}\t{smiles}\n")

            total_compounds += len(compounds)

        return total_compounds
