# chembl_curator/config.py

import json
from dataclasses import dataclass, asdict
from typing import List, Dict
from pathlib import Path


@dataclass
class CurationConfig:
    # Activity thresholds
    activity_thresholds: Dict[str, float] = None

    # Compound size thresholds
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

    # Validity filtering options
    require_standard_flag: bool = False  # Require curated/standardized data
    exclude_invalid_data: bool = True    # Exclude data with validity comments
    exclude_duplicates: bool = True      # Exclude potential duplicates

    # Assay quality filters
    min_confidence_score: int = None     # Minimum assay confidence (0-9, 9=highest)
    assay_types: List[str] = None        # Assay types (B=Binding, F=Functional, etc.)
    bao_formats: List[str] = None        # BAO format IDs (e.g., BAO_0000357)

    # Activity value filters
    min_pchembl_value: float = None      # Minimum pChEMBL value (-log10 molar activity)
    
    def __post_init__(self):
        """Set default values"""
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

        # Defaults for new filters (None means no filtering)
        # User can set these explicitly for stricter filtering
    
    @classmethod
    def from_file(cls, config_path: Path) -> 'CurationConfig':
        with open(config_path, 'r') as f:
            data = json.load(f)
        return cls(**data)
    
    def to_file(self, config_path: Path) -> None:
        with open(config_path, 'w') as f:
            json.dump(asdict(self), f, indent=2)
    
    @classmethod
    def create_example_config(cls, output_path: Path) -> None:
        example_config = cls()
        example_config.to_file(output_path)
