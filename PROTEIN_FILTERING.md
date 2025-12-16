# Protein Filtering Pipeline

This document describes the protein filtering pipeline that filters targets based on PDB structure availability and binding site analysis.

## Overview

The protein filtering pipeline processes curated ChEMBL targets and filters them based on:

1. **PDB Structure Availability**: Fetches PDB structures from UniProt
2. **AlphaFold Model Availability**: Downloads AlphaFold predicted structures
3. **Ligand-Bound Structures**: Identifies structures with biologically relevant ligands
4. **Single Binding Site**: Ensures all ligands across structures are in the same pocket

## Pipeline Steps

### Step 1: Fetch PDB Information
- Queries UniProt REST API for each target
- Retrieves PDB IDs, experimental method, resolution, and chain information
- Saves results to `pdbid.list` file in each target directory

### Step 2: Download PDB Structures
- Downloads PDB files using `pdb_get` command
- Creates `pdb/` directory under each target
- Saves structures in standard PDB format

### Step 3: Download AlphaFold Models
- Downloads AlphaFold predicted structure from EBI
- URL format: `https://alphafold.ebi.ac.uk/files/AF-{uniid}-F1-model_v6.pdb`
- Saves as `AF-{uniid}.pdb` in the pdb directory

### Step 4: Detect Ligand-Bound Structures
- Parses HETATM records from PDB files
- Excludes non-biological ligands (water, ions, buffers, crystallization aids)
- Checks if ligand has contact with protein (within 4Å of any heavy atom)
- Identifies single-ligand or clustered-ligand structures (ligands within 10Å)
- Saves list to `ligand_bound_pdbs.txt`

### Step 5: Align Structures
- Aligns all ligand-bound structures to AlphaFold model using TMalign
- Only saves the chain of interest from each structure
- Saves aligned structures to `aligned/` directory

### Step 6: Check Pocket Clustering
- Re-extracts ligand centers from aligned structures
- Checks if all ligands are in the same pocket (centers within 10Å)
- Saves pocket information to `pocket_info.csv` with columns:
  - PDB_ID
  - Ligand_Name
  - Center_X, Center_Y, Center_Z

### Step 7: Filter and Save Results
- Keeps only targets with single binding sites
- Saves passed target IDs to `passed_targets.txt`
- Optionally deletes filtered-out target directories

## Directory Structure

After running the protein filtering pipeline, each target directory will contain:

```
{UNIPROT_ID}/
├── comps/                    # From ligand filtering
│   └── smiles/
│       └── *.smi
├── pdb/                      # Downloaded structures
│   ├── AF-{UNIPROT_ID}.pdb  # AlphaFold model
│   ├── {pdbid}.pdb          # Experimental structures
│   └── ...
├── aligned/                  # Aligned structures
│   ├── {pdbid}_aligned.pdb
│   └── ...
├── pdbid.list               # PDB information
├── ligand_bound_pdbs.txt    # Ligand-bound PDB IDs
└── pocket_info.csv          # Pocket clustering info
```

## Usage

### Using CLI

```bash
# Run protein filtering on curated data
chembl-curator filter-proteins --curated-dir curated_data_filtered --n-processes 8

# With custom log level
chembl-curator filter-proteins --curated-dir curated_data_filtered --log-level DEBUG
```

### Using Python API

```python
from pathlib import Path
from chembl_curator import ProteinFilter

# Create protein filter
protein_filter = ProteinFilter(
    curated_dir=Path("curated_data_filtered"),
    log_level="INFO"
)

# Run pipeline (sequential)
passed_targets = protein_filter.run_pipeline(n_processes=1)

# Run pipeline (parallel with 8 processes)
passed_targets = protein_filter.run_pipeline(n_processes=8)

print(f"Passed {len(passed_targets)} targets")
```

### Using Example Script

```bash
python example_scripts/run_protein_filter.py \
    --curated-dir curated_data_filtered \
    --n-processes 8 \
    --log-level INFO
```

## Requirements

### External Tools
- `pdb_get`: PDB download tool (must be in PATH)
- `TMalign`: Structure alignment tool (expected at `/applic/bin/TMalign`)
- `wget`: For downloading AlphaFold models

### Python Dependencies
- requests: For UniProt API calls
- numpy: For coordinate calculations

## Error Handling

The pipeline is designed to be robust:

- **No AlphaFold model**: Target is filtered out, pipeline continues
- **No PDB structures**: Target is filtered out, pipeline continues
- **Download failures**: Logged and skipped
- **Alignment failures**: Structure is skipped, others continue
- **Parse errors**: Logged and structure is skipped

All errors are logged but don't stop the pipeline.

## Filtering Criteria

### Excluded Ligands

The following are excluded as non-biological ligands:
- **Sugars and glycans**: NAG, MAN, GLC, etc.
- **Ions**: CA, MG, ZN, NA, CL, etc.
- **Buffers**: MES, HEPES, TRIS, etc.
- **Crystallization aids**: PEG, GOL, EDO, SO4, etc.
- **Water**: HOH, H2O, WAT

See `exclusion.py` for full list.

### Contact Definition

A ligand is considered in contact with protein if:
- Any ligand heavy atom is within 4Å of any protein heavy atom
- Only considers heavy atoms (excludes hydrogens)

### Single Ligand Definition

A structure is considered single-ligand if either:
1. Only one ligand is present, OR
2. Multiple ligands are clustered (all centers within 10Å of each other)

### Same Pocket Definition

Ligands from different structures are in the same pocket if:
- All ligand centers are within 10Å of each other (after alignment)

## Output Files

### passed_targets.txt
List of UniProt IDs that passed all filters:
```
P28222
Q9Y6K9
P00533
...
```

### pdbid.list
PDB information for each target:
```
# PDBID method resolution chains
1ATP  X-RAY DIFFRACTION  2.0 A  A/B=1-250
2ATP  X-RAY DIFFRACTION  1.8 A  A=1-250
...
```

### ligand_bound_pdbs.txt
Ligand-bound PDB IDs:
```
1ATP  ATP
2ATP  ADP,MG
...
```

### pocket_info.csv
Pocket clustering information:
```
PDB_ID,Ligand_Name,Center_X,Center_Y,Center_Z
1ATP,ATP,12.345,23.456,34.567
2ATP,ADP,12.123,23.234,34.345
...
```

## Performance

- **Sequential processing**: ~30-60 seconds per target (depending on number of PDBs)
- **Parallel processing**: Linear speedup with number of cores
- **Recommended**: 8-16 processes for optimal performance

## Troubleshooting

### pdb_get not found
Ensure `pdb_get` is installed and in your PATH:
```bash
which pdb_get
```

### TMalign not found
Update the TMalign path in `protein_filter.py` or create symlink:
```bash
ln -s /path/to/TMalign /applic/bin/TMalign
```

### wget not available
Install wget:
```bash
# Ubuntu/Debian
sudo apt-get install wget

# macOS
brew install wget
```

### Memory issues with parallel processing
Reduce the number of processes:
```bash
chembl-curator filter-proteins --curated-dir curated_data_filtered --n-processes 4
```

## Examples

### Example 1: Test on single target
```python
from pathlib import Path
from chembl_curator.protein_filter import ProteinFilter

pf = ProteinFilter(Path("curated_data_filtered"), log_level="DEBUG")
result = pf.process_target("P28222")
print(f"Target P28222 passed: {result}")
```

### Example 2: Process specific targets
```python
from pathlib import Path
from chembl_curator.protein_filter import ProteinFilter

pf = ProteinFilter(Path("curated_data_filtered"))

targets = ["P28222", "Q9Y6K9", "P00533"]
passed = []

for target in targets:
    if pf.process_target(target):
        passed.append(target)

print(f"Passed: {passed}")
```

### Example 3: Custom distance thresholds
You can modify the distance thresholds in the code:
- `is_single_ligand_bound(ligands, distance_threshold=10.0)`: Ligand clustering
- `check_same_pocket(centers, distance_threshold=10.0)`: Pocket clustering

## Citation

If you use this pipeline, please cite:
- UniProt: https://www.uniprot.org
- AlphaFold: https://alphafold.ebi.ac.uk
- TMalign: Zhang & Skolnick (2005) Nucleic Acids Res. 33, 2302-9
- ChEMBL: https://www.ebi.ac.uk/chembl/
