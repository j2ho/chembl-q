
# ============================================================================
# chembl_curator/filters.py
# ============================================================================

"""Filtering utilities for ChEMBL data."""

from typing import Optional
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging

from .config import CurationConfig


class ActivityFilter:
    """Filter bioactivity data based on potency thresholds."""
    
    def __init__(self, config: CurationConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def is_active(self, value: float, units: str) -> bool:
        """Check if compound is active based on value and units.
        
        Args:
            value: Activity value
            units: Units (nM, uM, etc.)
            
        Returns:
            True if compound is considered active
        """
        threshold = self.config.activity_thresholds.get(units)
        if threshold is None:
            return False
        
        return value < threshold


class CompoundFilter:
    """Filter compounds based on chemical properties."""
    
    def __init__(self, config: CurationConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def passes_filters(self, smiles: str) -> bool:
        """Check if compound passes all filters.
        
        Args:
            smiles: SMILES string
            
        Returns:
            True if compound passes all filters
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            
            # Heavy atom count filter
            heavy_atoms = mol.GetNumHeavyAtoms()
            if not (self.config.min_heavy_atoms <= heavy_atoms <= self.config.max_heavy_atoms):
                return False
            
            return True
            
        except Exception as e:
            self.logger.warning(f"Error processing SMILES {smiles}: {e}")
            return False
    
    def get_properties(self, smiles: str) -> Optional[dict]:
        """Get molecular properties for a compound.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary of molecular properties or None if invalid
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            return {
                'heavy_atoms': mol.GetNumHeavyAtoms(),
                'mw': Descriptors.MolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'hbd': Descriptors.NumHDonors(mol),
                'hba': Descriptors.NumHAcceptors(mol),
                'tpsa': Descriptors.TPSA(mol)
            }
            
        except Exception as e:
            self.logger.warning(f"Error calculating properties for {smiles}: {e}")
            return None
