
# ============================================================================
# chembl_curator/curator.py
# ============================================================================

"""Main curation pipeline for ChEMBL data."""

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
    """Results from curation pipeline."""
    total_activities: int
    filtered_activities: int
    total_proteins: int
    total_compounds: int
    output_directory: Path


class ChEMBLCurator:
    """Main class for ChEMBL bioactivity data curation."""
    
    def __init__(
        self, 
        config: Optional[CurationConfig] = None,
        log_level: str = "INFO"
    ):
        """Initialize ChEMBL Curator.
        
        Args:
            config: Configuration for curation parameters
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
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
        """Download ChEMBL SQLite database.
        
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
        """Run complete curation pipeline.
        
        Args:
            database_path: Path to ChEMBL SQLite database
            output_dir: Output directory for curated data
            
        Returns:
            CurationResults with statistics
        """
        self.logger.info("Starting ChEMBL curation pipeline")
        
        # Download database if not provided
        if database_path is None:
            database_path = self.download_database()
            
        # Connect to database
        conn = sqlite3.connect(database_path)
        
        try:
            # Extract bioactivity data
            activities = self._extract_activities(conn)
            self.logger.info(f"Extracted {len(activities)} activity records")
            
            # Get mappings
            chembl_to_uniprot = self._get_uniprot_mapping(conn)
            chembl_to_smiles = self._get_smiles_mapping(conn)
            
            # Filter and organize
            results = self._process_activities(
                activities, 
                chembl_to_uniprot, 
                chembl_to_smiles, 
                output_dir
            )
            
            self.logger.info(f"Curation completed: {results.total_compounds} compounds for {results.total_proteins} proteins")
            return results
            
        finally:
            conn.close()
    
    def _extract_activities(self, conn: sqlite3.Connection) -> List[Tuple]:
        """Extract bioactivity data from database."""
        query = """
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
            td.target_type = 'SINGLE PROTEIN'
            AND a.standard_type IN ('Kd', 'Ki', 'IC50', 'EC50')
            AND a.standard_relation IN ('=', '<=')
            AND a.standard_units IN ('nM', 'uM')
            AND a.standard_value IS NOT NULL
        """
        return conn.execute(query).fetchall()
    
    def _get_uniprot_mapping(self, conn: sqlite3.Connection) -> Dict[str, str]:
        """Get ChEMBL to UniProt mapping."""
        query = """
        SELECT DISTINCT 
            td.chembl_id AS target_chembl_id,
            cs.accession AS uniprot_id
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE 
            td.target_type = 'SINGLE PROTEIN'
            AND cs.accession IS NOT NULL
            AND (cs.accession LIKE 'P%' OR cs.accession LIKE 'Q%' OR cs.accession LIKE 'O%')
        """
        results = conn.execute(query).fetchall()
        return {chembl_id: uniprot_id for chembl_id, uniprot_id in results}
    
    def _get_smiles_mapping(self, conn: sqlite3.Connection) -> Dict[str, str]:
        """Get compound SMILES mapping."""
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
        """Process and filter activities."""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        protein_compounds = {}
        filtered_count = 0
        
        for activity in activities:
            target_chembl, compound_chembl, std_type, std_relation, std_value, std_units = activity
            
            # Apply activity filter
            if not self.activity_filter.is_active(std_value, std_units):
                continue
                
            # Map to UniProt
            if target_chembl not in chembl_to_uniprot:
                continue
            uniprot_id = chembl_to_uniprot[target_chembl]
            
            # Get SMILES
            if compound_chembl not in chembl_to_smiles:
                continue
            smiles = chembl_to_smiles[compound_chembl]
            
            # Apply compound filters
            if not self.compound_filter.passes_filters(smiles):
                continue
            
            # Store compound
            if uniprot_id not in protein_compounds:
                protein_compounds[uniprot_id] = set()
            protein_compounds[uniprot_id].add((compound_chembl, smiles))
            filtered_count += 1
        
        # Write output files
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
        """Write compound files organized by protein."""
        total_compounds = 0
        
        for uniprot_id, compounds in protein_compounds.items():
            prot_dir = output_dir / uniprot_id / "comps"
            prot_dir.mkdir(parents=True, exist_ok=True)
            
            for compound_id, smiles in compounds:
                smi_file = prot_dir / f"{compound_id}.smi"
                with open(smi_file, 'w') as f:
                    f.write(smiles)
                total_compounds += 1
        
        return total_compounds
