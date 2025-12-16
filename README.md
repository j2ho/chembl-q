# ChEMBL Curator

A comprehensive pipeline for curating ChEMBL database and filtering protein-ligand complexes based on activity data, compound properties, and structural quality.

## Overview

ChEMBL Curator provides a two-stage pipeline:

1. **Compound Filtering Stage**: Extracts and filters ligands from ChEMBL database based on:
   - Activity data (Ki, Kd, IC50, EC50)
   - Activity thresholds and relations
   - Compound validity (SMILES, molecular properties)
   - Target types (single proteins)

2. **Protein Filtering Stage**: Filters targets based on structural data:
   - PDB structure availability
   - AlphaFold model availability
   - Ligand-bound structures
   - Single binding site requirement

## Features

- **Flexible Configuration**: JSON-based configuration or command-line options
- **Robust Filtering**: Multiple validity checks for activity data and compounds
- **Structure-Based Filtering**: Automated PDB download, alignment, and pocket analysis
- **Parallel Processing**: Support for multi-core processing
- **Comprehensive Logging**: Detailed logging at multiple levels
- **Error Handling**: Continues processing on errors, logs issues

## Installation

### Prerequisites

- Python 3.8+
- External tools (for protein filtering):
  - `pdb_get`: PDB download tool
  - `TMalign`: Structure alignment tool (expected at `./bin/TMalign`)
  - `wget`: For downloading AlphaFold models

### Python Dependencies

```bash
# Install the package
pip install -e .

# Or install dependencies manually
pip install click requests numpy rdkit pandas
```

### External Tools Setup

```bash
# Check if tools are available
which pdb_get
which wget
ls ./bin/TMalign

# If TMalign is in a different location, update protein_filter.py or create symlink
mkdir -p ./bin
ln -s /path/to/TMalign ./bin/TMalign
```

## Quick Start

### Complete Pipeline (Compound + Protein Filtering)

```bash
# Step 1: Compound filtering - creates curated_data_filtered/
echo "Step 1: Running compound filtering..."
chembl-curator curate --download --output curated_data_filtered --log-level INFO

# Wait for completion, then:

# Step 2: Protein filtering - adds PDB structures and filters by binding site
echo "Step 2: Running protein filtering..."
chembl-curator filter-proteins \
    --curated-dir curated_data_filtered \
    --n-processes 8 \
    --log-level INFO

# Step 3: Check results
echo "Step 3: Checking results..."
cat curated_data_filtered/passed_targets.txt
echo "Total passed targets: $(wc -l < curated_data_filtered/passed_targets.txt)"
```

### Testing on Subset (Recommended First)

For testing, you might want to try on a small subset first:

```bash
# Step 1: Create test directory with just a few targets
mkdir -p test_curated
cp -r curated_data_filtered/P28222 test_curated/
cp -r curated_data_filtered/Q9Y6K9 test_curated/
cp -r curated_data_filtered/P00533 test_curated/

# Step 2: Run protein filtering on test subset
chembl-curator filter-proteins \
    --curated-dir test_curated \
    --n-processes 1 \
    --log-level DEBUG

# Step 3: Check if it worked
ls -la test_curated/P28222/pdb/
cat test_curated/passed_targets.txt
```

## Usage

### Stage 1: Compound Filtering

#### Using CLI

```bash
# Basic usage - download and curate with defaults
chembl-curator curate --download --output curated_data_filtered

# Using existing database
chembl-curator curate --database /path/to/chembl_34.db --output curated_data_filtered

# Using custom configuration file
chembl-curator curate --download --config config.json --output curated_data_filtered

# Override specific parameters
chembl-curator curate --download --output curated_data_filtered \
    --activity-types Ki Kd IC50 \
    --relations "=" "<=" "<" \
    --units nM

# Create example config file
chembl-curator curate --create-config config.json

# With debug logging
chembl-curator curate --download --output curated_data_filtered --log-level DEBUG
```

#### Using Python API

```python
from chembl_curator import ChEMBLCurator, CurationConfig
from pathlib import Path

# Option 1: Using default configuration
curator = ChEMBLCurator(log_level='INFO')
results = curator.run_pipeline(
    database_path=None,  # Will download ChEMBL
    output_dir=Path('curated_data_filtered')
)

print(f"Total compounds: {results.total_compounds}")
print(f"Total proteins: {results.total_proteins}")

# Option 2: Using custom configuration
config = CurationConfig(
    activity_thresholds={'nM': 1000.0, 'uM': 1.0},
    activity_types=['Ki', 'Kd', 'IC50'],
    units=['nM'],
    min_pchembl_value=6.0,
    min_confidence_score=8
)

curator = ChEMBLCurator(config=config, log_level='INFO')
results = curator.run_pipeline(
    database_path=Path('chembl_34.db'),
    output_dir=Path('curated_data_filtered')
)
```

### Stage 2: Protein Filtering

#### Using CLI

```bash
# Sequential processing (1 core)
chembl-curator filter-proteins --curated-dir curated_data_filtered --n-processes 1

# Parallel processing (8 cores, recommended)
chembl-curator filter-proteins --curated-dir curated_data_filtered --n-processes 8

# With debug logging
chembl-curator filter-proteins --curated-dir curated_data_filtered --n-processes 8 --log-level DEBUG
```

#### Using Python API

```python
from chembl_curator import ProteinFilter
from pathlib import Path

# Create protein filter
pf = ProteinFilter(
    curated_dir=Path('curated_data_filtered'),
    log_level='INFO'
)

# Run pipeline with parallel processing
passed_targets = pf.run_pipeline(n_processes=8)

print(f"Passed {len(passed_targets)} targets with single binding sites")

# Process single target
result = pf.process_target('P28222')
print(f"Target P28222 passed: {result}")
```

#### Using Example Script

```bash
python example_scripts/run_protein_filter.py \
    --curated-dir curated_data_filtered \
    --n-processes 8 \
    --log-level INFO
```

## Configuration

### Default Configuration

The default configuration filters for:

**Activity Types:**
- `Ki`, `Kd`, `IC50`, `EC50`

**Relations:**
- `=`, `<=`

**Units:**
- `nM` (nanomolar), `uM` (micromolar)

**Activity Thresholds:**
- ≤ 10000 nM (10 µM) for nM units
- ≤ 10 µM for µM units

**Target Types:**
- `SINGLE PROTEIN` only

**Compound Filters:**
- Heavy atoms: 5-80
- Valid SMILES required

**Validity Filters:**
- Standard flag: not required (false)
- Exclude invalid data: true
- Exclude duplicates: true
- Minimum confidence score: ≥ 6
- Minimum pChEMBL value: ≥ 5.0

**Assay Types:**
- `B` (Binding assays only)

**BAO Formats:**
- `BAO_0000357` (single protein format)

### Custom Configuration File

Create a JSON configuration file:

```json
{
  "activity_thresholds": {
    "nM": 1000.0,
    "uM": 1.0
  },
  "min_heavy_atoms": 5,
  "max_heavy_atoms": 50,
  "target_types": ["SINGLE PROTEIN"],
  "activity_types": ["Ki", "Kd", "IC50"],
  "relations": ["=", "<="],
  "units": ["nM", "uM"],
  "require_standard_flag": false,
  "exclude_invalid_data": true,
  "exclude_duplicates": true,
  "min_confidence_score": 8,
  "assay_types": ["B", "F"],
  "bao_formats": ["BAO_0000357"],
  "min_pchembl_value": 6.0
}
```

Use with:
```bash
chembl-curator curate --config config.json --download --output curated_data_filtered
```

## Output Structure

### After Compound Filtering

```
curated_data_filtered/
├── P28222/                      # UniProt ID
│   └── comps/
│       └── smiles/
│           ├── CHEMBL31217.smi
│           ├── CHEMBL308727.smi
│           └── ...              # One .smi file per compound
├── Q9Y6K9/
│   └── comps/
│       └── smiles/
│           └── ...
└── ...                          # One directory per target
```

### After Protein Filtering

```
curated_data_filtered/
├── P28222/                      # Target that passed all filters
│   ├── comps/
│   │   └── smiles/
│   │       └── *.smi
│   ├── pdb/                     # Downloaded structures
│   │   ├── AF-P28222.pdb       # AlphaFold model
│   │   ├── 1atp.pdb            # Experimental PDB structures
│   │   ├── 2atp.pdb
│   │   └── ...
│   ├── aligned/                 # Structures aligned to AF model (per chain)
│   │   ├── 1atp_A.pdb
│   │   ├── 1atp_B.pdb
│   │   ├── 2atp_A.pdb
│   │   └── ...
│   ├── pdbid.list              # PDB information
│   ├── ligand_bound_pdbs.txt   # Ligand-bound structures
│   └── pocket_info.csv         # Pocket clustering data
├── Q9Y6K9/                      # Another target that passed
│   └── ...
├── passed_targets.txt           # List of targets that passed protein filtering
└── ...
```

### Output Files Description

**pdbid.list** - PDB structure information:
```
# PDBID method resolution chains
1ATP  X-RAY DIFFRACTION  2.0 A  A/B=1-250
2ATP  X-RAY DIFFRACTION  1.8 A  A=1-250
```

**ligand_bound_pdbs.txt** - Structures with biological ligands:
```
1ATP  Chain A: ATP
1ATP  Chain B: ATP
2ATP  Chain A: ADP,MG
```

**pocket_info.csv** - Ligand pocket coordinates (from aligned structures):
```csv
PDB_ID,Chain,Aligned_File,Ligand_Name,Center_X,Center_Y,Center_Z
1ATP,A,1atp_A.pdb,ATP,12.345,23.456,34.567
1ATP,B,1atp_B.pdb,ATP,12.123,23.234,34.345
2ATP,A,2atp_A.pdb,ADP,12.100,23.200,34.300
```

**passed_targets.txt** - Final filtered targets:
```
P28222
Q9Y6K9
P00533
```

## Protein Filtering Pipeline Details

The protein filtering stage consists of 7 steps:

### 1. Fetch PDB Information
- Queries UniProt REST API
- Retrieves PDB ID, method, resolution, chain information

### 2. Download PDB Structures
- Uses `pdb_get` command
- Creates `pdb/` directory

### 3. Download AlphaFold Models
- Downloads from `https://alphafold.ebi.ac.uk/files/AF-{uniid}-F1-model_v6.pdb`
- Used as reference for alignment

### 4. Detect Ligand-Bound Structures
- Parses HETATM records
- Excludes non-biological ligands (ions, buffers, crystallization aids, water)
- Checks ligand-protein contact (≤4Å)
- Identifies single-ligand or clustered structures (≤10Å between ligand centers)

### 5. Align Structures
- Aligns to AlphaFold model using TMalign
- Saves only chain of interest

### 6. Check Pocket Clustering
- Re-extracts ligand centers from aligned structures
- Verifies all ligands in same pocket (≤10Å)

### 7. Filter and Save
- Keeps only single binding site targets
- Saves results to `passed_targets.txt`

### Excluded Ligands

Non-biological ligands that are filtered out:

- **Water**: HOH, H2O, WAT
- **Ions**: CA, MG, ZN, NA, CL, FE, etc.
- **Buffers**: MES, HEPES, TRIS, etc.
- **Sugars**: NAG, MAN, GLC, GAL, etc.
- **Crystallization aids**: PEG, GOL, EDO, SO4, etc.

See `example_scripts/exclusion.py` or `chembl_curator/protein_filter.py` for complete list.

## Examples

### Example 1: Basic Usage

```python
from chembl_curator import ChEMBLCurator, ProteinFilter
from pathlib import Path

# Stage 1: Compound filtering
curator = ChEMBLCurator(log_level='INFO')
results = curator.run_pipeline(output_dir=Path('curated_data'))

# Stage 2: Protein filtering
pf = ProteinFilter(curated_dir=Path('curated_data'), log_level='INFO')
passed = pf.run_pipeline(n_processes=8)

print(f"Final database: {len(passed)} targets")
```

### Example 2: Custom Thresholds

```python
from chembl_curator import ChEMBLCurator, CurationConfig
from pathlib import Path

# More stringent activity threshold
config = CurationConfig(
    activity_thresholds={'nM': 100.0, 'uM': 0.1},  # 100 nM instead of 10 µM
    activity_types=['Ki', 'Kd'],
    units=['nM', 'uM'],
    min_pchembl_value=7.0,      # Higher quality data
    min_confidence_score=8,      # Higher confidence
    max_heavy_atoms=30,          # Smaller molecules
    assay_types=['B']            # Binding assays only
)

curator = ChEMBLCurator(config=config, log_level='INFO')
results = curator.run_pipeline(output_dir=Path('curated_high_quality'))
```

### Example 3: Process Specific Targets

```python
from chembl_curator.protein_filter import ProteinFilter
from pathlib import Path

pf = ProteinFilter(Path('curated_data_filtered'))

# Process specific targets
targets = ['P28222', 'Q9Y6K9', 'P00533']
passed = []

for target in targets:
    try:
        if pf.process_target(target):
            passed.append(target)
            print(f"✓ {target} passed")
        else:
            print(f"✗ {target} filtered out")
    except Exception as e:
        print(f"✗ {target} error: {e}")

print(f"\nPassed: {passed}")
```

### Example 4: Inspect Results

```python
from pathlib import Path
import pandas as pd

# Read passed targets
with open('curated_data_filtered/passed_targets.txt') as f:
    targets = [line.strip() for line in f]

print(f"Total passed targets: {len(targets)}")

# Inspect a specific target
target = targets[0]
target_dir = Path('curated_data_filtered') / target

# Count compounds
compounds = list((target_dir / 'comps' / 'smiles').glob('*.smi'))
print(f"\n{target}: {len(compounds)} compounds")

# Read PDB list
with open(target_dir / 'pdbid.list') as f:
    pdbs = [line for line in f if not line.startswith('#')]
print(f"PDB structures: {len(pdbs)}")

# Read pocket info
pocket_df = pd.read_csv(target_dir / 'pocket_info.csv')
print(f"\nPocket info:")
print(pocket_df)
```

## Performance

### Compound Filtering
- **Time**: 10-30 minutes (depending on database size)
- **Memory**: ~2-4 GB
- **Output**: ~1000-5000 targets (depends on filters)

### Protein Filtering
- **Per target**: 30-60 seconds
- **Total (4000 targets, 8 cores)**: 3-6 hours
- **Memory**: 100-500 MB per process
- **Recommended**: 8-16 processes for optimal performance

### Disk Space

- **ChEMBL database**: ~3 GB
- **Compound filtering output**: ~10-50 MB
- **Protein filtering output**: ~500 MB - 5 GB (depends on number of PDBs)

## Troubleshooting

### Common Issues

**1. pdb_get not found**
```bash
# Check if installed
which pdb_get

# If not, install or add to PATH
export PATH=/path/to/pdb_tools:$PATH
```

**2. TMalign not found**
```bash
# Create bin directory and symlink to TMalign
mkdir -p ./bin
ln -s /path/to/TMalign ./bin/TMalign

# Or update path in protein_filter.py line 460
```

**3. Memory issues with parallel processing**
```bash
# Reduce number of processes
chembl-curator filter-proteins --curated-dir curated_data_filtered --n-processes 4
```

**4. Download failures**
- Check internet connection
- Some PDBs may be obsolete (pipeline continues)
- Some targets may not have AlphaFold models (filtered out)

**5. Import errors**
```bash
# Reinstall package
pip install -e .

# Or check Python path
export PYTHONPATH=/path/to/ChEMBL-Q:$PYTHONPATH
```

### Debug Mode

For detailed error messages:

```bash
# Compound filtering
chembl-curator curate --download --output test_output --log-level DEBUG

# Protein filtering
chembl-curator filter-proteins --curated-dir test_output --n-processes 1 --log-level DEBUG
```

## Project Structure

```
ChEMBL-Q/
├── chembl_curator/              # Main package
│   ├── __init__.py
│   ├── cli.py                   # Command-line interface
│   ├── config.py                # Configuration classes
│   ├── curator.py               # Compound filtering
│   ├── downloader.py            # ChEMBL downloader
│   ├── filters.py               # Activity and compound filters
│   ├── protein_filter.py        # Protein filtering pipeline
│   └── utils.py                 # Utility functions
├── example_scripts/             # Example usage scripts
│   ├── getpdb.py
│   ├── exclusion.py
│   ├── aligntm.py
│   └── run_protein_filter.py
├── pyproject.toml               # Package configuration
├── README.md                    # This file
├── PROTEIN_FILTERING.md         # Detailed protein filtering docs
└── .gitignore
```

## Documentation

- **README.md** (this file): Overview and quick start guide
- **PROTEIN_FILTERING.md**: Detailed protein filtering pipeline documentation
- **example_scripts/**: Example usage scripts and reference implementations

## Citation

If you use ChEMBL Curator in your research, please cite:

- **ChEMBL**: Mendez D, et al. (2019) Nucleic Acids Res. 47(D1):D930-D940
- **UniProt**: The UniProt Consortium (2023) Nucleic Acids Res. 51(D1):D523-D531
- **AlphaFold**: Jumper J, et al. (2021) Nature 596:583-589
- **TMalign**: Zhang Y, Skolnick J (2005) Nucleic Acids Res. 33:2302-2309

## License

This project is provided as-is for research purposes.

## Author

Jiho Sim

## Acknowledgments

- ChEMBL database from EMBL-EBI
- AlphaFold structures from DeepMind/EMBL-EBI
- PDB structures from RCSB/PDBe
- RDKit for cheminformatics
- TMalign for structure alignment
