# ============================================================================
# README.md
# ============================================================================

"""
# ChEMBL Curator

A Python package for curating bioactivity data from the ChEMBL database.

## Features

- Automated ChEMBL database download and setup
- Configurable activity and compound filtering
- UniProt target mapping
- Organized output by protein target
- Command-line interface
- Comprehensive logging

## Installation

```bash
pip install chembl-curator
```

Or install from source:

```bash
git clone https://github.com/username/chembl-curator.git
cd chembl-curator
pip install -e .
```

## Quick Start

### Python API

```python
from chembl_curator import ChEMBLCurator

# Create curator with default settings
curator = ChEMBLCurator()

# Download database and run curation
results = curator.run_pipeline(output_dir="./my_curated_data")

print(f"Curated {results.total_compounds} compounds for {results.total_proteins} proteins")
```

### Command Line

```bash
# Download database and run curation
chembl-curator --download --output ./curated_data

# Use existing database
chembl-curator --database ./chembl_35.db --output ./curated_data
```

## Configuration

Customize curation parameters:

```python
from chembl_curator import ChEMBLCurator, CurationConfig

config = CurationConfig(
    activity_thresholds={'nM': 1000, 'uM': 1},  # Custom thresholds
    min_heavy_atoms=10,                          # Larger molecules only
    max_heavy_atoms=50,                          # Smaller molecules only
    activity_types=['IC50', 'Ki']                # Specific assay types
)

curator = ChEMBLCurator(config=config)
results = curator.run_pipeline()
```

## Output Structure

```
output_dir/
├── P12345/                 # UniProt ID
│   └── comps/
│       ├── CHEMBL123.smi   # SMILES files
│       ├── CHEMBL456.smi
│       └── ...
└── Q67890/                 # Another protein
    └── comps/
        └── ...
```

## Requirements

- Python >= 3.8
- RDKit
- pandas
- requests
- tqdm
- click

## License

MIT License

## Citation

If you use this package in your research, please cite:

```
@software{chembl_curator,
  title={ChEMBL Curator: A Python package for ChEMBL bioactivity data curation},
  author={Your Name},
  url={https://github.com/username/chembl-curator},
  year={2024}
}
```
"""

