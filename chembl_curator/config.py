
# ============================================================================
# chembl_curator/config.py
# ============================================================================

"""Configuration classes for ChEMBL curation."""

from dataclasses import dataclass
from typing import List, Dict


@dataclass
class CurationConfig:
    """Configuration for ChEMBL curation parameters."""
    
    # Activity thresholds
    activity_thresholds: Dict[str, float] = None
    
    # Compound filters
    min_heavy_atoms: int = 5
    max_heavy_atoms: int = 80
    
    # Target types to include
    target_types: List[str] = None
    
    # Activity types to include  
    activity_types: List[str] = None
    
    # Relations to include
    relations: List[str] = None
    
    # Units to include
    units: List[str] = None
    
    def __post_init__(self):
        """Set default values."""
        if self.activity_thresholds is None:
            self.activity_thresholds = {
                'nM': 10000.0,
                'uM': 10.0
            }
        
        if self.target_types is None:
            self.target_types = ['SINGLE PROTEIN']
            
        if self.activity_types is None:
            self.activity_types = ['Kd', 'Ki', 'IC50', 'EC50']
            
        if self.relations is None:
            self.relations = ['=', '<=']
            
        if self.units is None:
            self.units = ['nM', 'uM']
