# chembl_curator/curator.py

import sqlite3
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from .downloader import ChEMBLDownloader
from .filters import ActivityFilter, CompoundFilter
from .config import CurationConfig
from .utils import setup_logging


@dataclass
class CurationResults:
    total_activities: int
    filtered_activities: int
    total_proteins: int
    total_compounds: int
    output_directory: Path


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
            
            results = self._process_activities(
                activities, 
                chembl_to_uniprot, 
                chembl_to_smiles, 
                output_dir
            )
            
            self.logger.info(f"Curation completed: {results.total_compounds} actives for {results.total_proteins} targets")
            return results
            
        finally:
            conn.close()
    
    def _extract_activities(self, conn: sqlite3.Connection) -> List[Tuple]:
        target_types = ','.join(f"'{t}'" for t in self.config.target_types)
        activity_types = ','.join(f"'{t}'" for t in self.config.activity_types)
        relations = ','.join(f"'{r}'" for r in self.config.relations)
        units = ','.join(f"'{u}'" for u in self.config.units)

        # Build validity filter conditions
        validity_conditions = []

        if self.config.require_standard_flag:
            validity_conditions.append("a.standard_flag = 1")

        if self.config.exclude_invalid_data:
            validity_conditions.append("a.data_validity_comment IS NULL")

        if self.config.exclude_duplicates:
            validity_conditions.append("(a.potential_duplicate IS NULL OR a.potential_duplicate = 0)")

        # Add assay quality filters
        if self.config.min_confidence_score is not None:
            validity_conditions.append(f"ass.confidence_score >= {self.config.min_confidence_score}")

        if self.config.assay_types is not None:
            assay_type_list = ','.join(f"'{t}'" for t in self.config.assay_types)
            validity_conditions.append(f"ass.assay_type IN ({assay_type_list})")

        if self.config.bao_formats is not None:
            bao_format_list = ','.join(f"'{b}'" for b in self.config.bao_formats)
            validity_conditions.append(f"ass.bao_format IN ({bao_format_list})")

        # Add pChEMBL filter
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
            a.standard_units
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
    ) -> CurationResults:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        protein_compounds = {}
        filtered_count = 0
        
        for activity in activities:
            target_chembl, compound_chembl, std_type, std_relation, std_value, std_units = activity
            
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
            
            if uniprot_id not in protein_compounds:
                protein_compounds[uniprot_id] = set()
            protein_compounds[uniprot_id].add((compound_chembl, smiles))
            filtered_count += 1
        
        total_compounds = self._write_output_files(protein_compounds, output_dir)
        
        return CurationResults(
            total_activities=len(activities),
            filtered_activities=filtered_count,
            total_proteins=len(protein_compounds),
            total_compounds=total_compounds,
            output_directory=output_dir
        )
    
    def _write_output_files(
        self, 
        protein_compounds: Dict[str, set], 
        output_dir: Path
    ) -> int:
        total_compounds = 0
        
        for uniprot_id, compounds in protein_compounds.items():
            prot_dir = output_dir / uniprot_id / "comps" / "smiles"
            prot_dir.mkdir(parents=True, exist_ok=True)
            
            for compound_id, smiles in compounds:
                smi_file = prot_dir / f"{compound_id}.smi"
                with open(smi_file, 'w') as f:
                    f.write(smiles)
                total_compounds += 1
        
        return total_compounds
