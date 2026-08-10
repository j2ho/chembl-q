<div align="center">
  <div>&nbsp;</div>
  <img src="docs/chemblq.png" width="300"/>

  <a href="https://j2ho.github.io/chembl-q/">Project Webpage</a>
</div>

A pipeline for curating ChEMBL into a virtual screening dataset for deep-learning model training/validation.

We filter ChEMBL activity data, select artificial decoys with reasonable compound/target criteria. Our pipeline also provides a leakage-resistant train/test splitting by receptor similarity. Using PDBbind+BioLip as default training set, easily replacable with user-custom bulk fasta. 

Using just our default setting, you can create the 'ChEMBL-LR' benchmark to test DL models trained on PDBbind and/or BioLip2 dataset wihtout having to worry about target leakage! 

---

## Pipeline Overview

| Stage | Command | Output |
|-------|---------|--------|
| 1. Compound filtering | `curate` | `{target}/actives.tsv`, `{target}/comps/smiles/*.smi` |
| 2. Protein filtering | `filter-proteins` | `aligned/*.pdb`, `pocket_info.csv`, `sequences.fasta`, `best_structure.tsv` |
| 3. Active clustering | `cluster-actives` | `{target}/actives_clustered.tsv` |
| 4. Compound pool | `build-pool` | `compound_pool.pkl` |
| 5. Receptor similarity | `receptor-sim` | `pairwise_seqid.tsv`, `pairwise_pocket_rmsd.tsv` |
| 6. Decoy selection | `select-decoys` | `{target}/decoys.tsv` |
| 7. Train/test split | `split` | `train.txt`, `test.txt` |

All outputs go under a single data directory (e.g. `curated_data_filtered/`).

---

## Installation

### Conda environment (recommended)

conda is the recommended approach - RDKit and nurikit are C++ extension packages that conda resolves cleanly.

```bash
conda create -n chemblq python=3.10
conda activate chemblq

# RDKit via conda-forge (simpler than pip for RDKit)
conda install -c conda-forge rdkit

# Remaining Python dependencies
pip install -e .          # installs click, numpy, pandas, requests, tqdm, nurikit
```

> **`nuri` vs `nurikit`**: the Python import is `import nuri` but the package name on PyPI is `nurikit`. The `pip install -e .` above handles this via the dependency in `pyproject.toml`.

### External binaries

These are not pip/conda packages and must be installed separately:

| Tool | Required for | Install |
|------|-------------|---------|
| **MMseqs2** | Stages 5, 7 (sequence search/clustering) | [github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2) - must be in `PATH` |
| `wget` | Stage 2 (AlphaFold download) | usually pre-installed |
| `pdb_get` | Stage 2 (optional) | local PDB mirror; falls back to RCSB web download |

> **Note:** Structure alignment (Stages 2, 5) uses `nurikit` (Python TMAlign bindings), which is installed automatically via `pip install -e .`. No separate TMalign binary is needed.

```bash
# Verify MMseqs2 is in PATH
mmseqs --help
```

---

## Quick Start

```bash
DATA=curated_data_filtered

# Stage 1: curate compounds from ChEMBL (bundled config.json applies opt-in filters)
chembl-curator curate --download --config config.json --output $DATA

# Stage 2: validate protein structures and binding sites
chembl-curator filter-proteins --curated-dir $DATA --n-processes 8

# Stage 3: cluster actives per target (Butina, Tanimoto ≥ 0.7)
chembl-curator cluster-actives --data-dir $DATA --workers 8

# Stage 4: build global compound pool
chembl-curator build-pool --data-dir $DATA

# Stage 5: compute pairwise receptor similarity
chembl-curator receptor-sim --data-dir $DATA --mode both --workers 8

# Stage 6: select property-matched, receptor-aware decoys
chembl-curator select-decoys --data-dir $DATA --max-decoys 30

# Stage 7: train/test split by sequence identity clustering
chembl-curator split --data-dir $DATA --valid-frac 0.1
```

---

## Stage Reference

### Stage 1: `curate`

Extracts and filters ligand-target pairs from ChEMBL.

```bash
chembl-curator curate --download --output curated_data_filtered
chembl-curator curate --database /path/to/chembl.db --config config.json --output curated_data_filtered
chembl-curator curate --create-config config.json   # generate example config
```

**Always on (built-in defaults):**
- Target type: SINGLE PROTEIN
- Activity types: Ki, Kd, IC50, EC50
- Relations: `=`, `<=`
- Units: nM, uM (< 10,000 nM / < 10 µM)
- Heavy atoms: 5-80, valid SMILES required
- Excludes rows with `data_validity_comment` or `potential_duplicate`

**Opt-in (only applied when `-c config.json` is passed):**
- `min_pchembl_value: 5.0`
- `min_confidence_score: 6`
- `assay_types: ["B"]` (binding only)
- `bao_formats: ["BAO_0000357"]`
- `require_standard_flag: true` (curated data only)

The bundled `config.json` at repo root enables all of the opt-in filters — pass it via `--config config.json` to reproduce the shipped dataset.

**Configuration (JSON):**
```json
{
  "activity_thresholds": {"nM": 10000.0, "uM": 10.0},
  "activity_types": ["Kd", "Ki", "IC50", "EC50"],
  "relations": ["=", "<="],
  "units": ["nM", "uM"],
  "min_pchembl_value": 5.0,
  "min_confidence_score": 6,
  "assay_types": ["B"]
}
```

---

### Stage 2: `filter-proteins`

Fetches PDB structures, downloads AlphaFold models, aligns structures, and keeps only targets with a single binding site.

```bash
chembl-curator filter-proteins --curated-dir curated_data_filtered --n-processes 8
```

Steps: fetch UniProt PDB list -> download PDB/AlphaFold -> detect ligand-bound structures -> align to AlphaFold (TMalign) -> cluster pockets -> filter single-site targets.

Also writes `sequences.fasta` (canonical UniProt sequences) and `best_structure.tsv` (best-resolution ligand-bound structure per target) needed by later stages.

---

### Stage 3: `cluster-actives`

Butina clustering of actives per target. Picks the highest-pChEMBL representative per cluster.

```bash
chembl-curator cluster-actives --data-dir curated_data_filtered --dist-thresh 0.3 --workers 8
```

| Option | Default | Description |
|--------|---------|-------------|
| `--dist-thresh` | 0.3 | Tanimoto distance threshold (0.3 -> similarity ≥ 0.7) |
| `--workers` | 1 | Parallel worker processes |

---

### Stage 4: `build-pool`

Builds a global compound pool from all clustered actives. Deduplicates by ChEMBL ID, computes molecular properties and Morgan fingerprints.

```bash
chembl-curator build-pool --data-dir curated_data_filtered
```

Pool contains: MW, cLogP, TPSA, HBD, HBA, aromatic rings, 2048-bit Morgan FP (radius 2), target membership set.

---

### Stage 5: `receptor-sim`

Computes pairwise receptor similarity via two independent methods.

```bash
# Both (recommended)
chembl-curator receptor-sim --data-dir curated_data_filtered --mode both --workers 8

# Sequence identity only (faster, no nuri required)
chembl-curator receptor-sim --data-dir curated_data_filtered --mode seqid
```

| Option | Default | Description |
|--------|---------|-------------|
| `--mode` | both | `seqid`, `pocket`, or `both` |
| `--seqid-threads` | 4 | MMseqs2 thread count |
| `--workers` | 4 | Processes for pocket RMSD |
| `--pocket-radius` | 10.0 | Pocket radius in Å |

**Sequence identity:** MMseqs2 all-vs-all -> `pairwise_seqid.tsv` (query, target, seqid 0-1)

**Pocket RMSD:** TM-align full chain -> filter to pocket residues within radius -> RMSD on ≥3 matched pairs -> `pairwise_pocket_rmsd.tsv`

---

### Stage 6: `select-decoys`

Selects property-matched, chemically dissimilar decoys per active. Excludes compounds active against receptors similar to the query target.

```bash
chembl-curator select-decoys --data-dir curated_data_filtered --max-decoys 30
```

| Option | Default | Description |
|--------|---------|-------------|
| `--max-decoys` | 30 | Decoys per active |
| `--seqid-thresh` | 0.6 | Seqid threshold for receptor exclusion |
| `--pocket-rmsd-thresh` | 3.0 | Pocket RMSD threshold (Å) for receptor exclusion |
| `--exclusion-mode` | or | `or` = exclude if seqid OR pocket matches |
| `--tanimoto-thresh` | 0.3 | Max Tanimoto between active and decoy |

Property matching windows: ±50 Da MW, ±2 cLogP, ±50 Å² TPSA, ±2 HBD, ±2 HBA, ±1 aromatic ring.

---

### Stage 7: `split`

Train/test split by sequence-identity clustering (MMseqs2). Greedy assignment balances per-source ratios. Per-entry sampling weight = `1 / log2(cluster_size + 1)`.

A bundled PDBbind+BioLip FASTA (`chembl_curator/assets/external_targets.fasta`) is included by default.

```bash
# Default: includes bundled PDBbind+BioLip sequences
chembl-curator split --data-dir curated_data_filtered --valid-frac 0.1

# ChEMBL-only split (no external sequences)
chembl-curator split --data-dir curated_data_filtered --no-external

# Use your own external FASTA
chembl-curator split --data-dir curated_data_filtered \
    --external-fasta /path/to/your.fasta
```

**External FASTA ID format:** IDs must be dot-prefixed as `>{source}.{entry_id}`, where `source` is any label (e.g. `pdbbind`, `biolip`, `myscreendb`) and `entry_id` is any identifier without spaces. For example:
```
>pdbbind.1a4k
MPPYTVVY...
>biolip.10gs_VWW_A_1
PYTVVYFP...
>myscreendb.custom_entry_42
MKWVTFIS...
```

**Output: `train.txt` / `test.txt`** (tab-separated, with header):
```
source	entry_id	compound	weight
chembl	A0A0H2UPP7	CHEMBL405346	0.17
chembl	A0A0H2UPP7	CHEMBL407216	0.17
biolip	10gs_VWW_A_1	-	0.19
pdbbind	1a4k	-	0.26
```

**Output: `chembl_targets.tsv`** (tab-separated, with header):
```
uniprot	split	n_actives	n_decoys
A0A0H2UPP7	train	2	60
Q9NR56	test	9	270
```

---

## Output Structure

```
curated_data_filtered/
├── sequences.fasta             # Canonical sequences, all passed targets
├── best_structure.tsv          # uniprot -> best PDB chain + resolution
├── compound_pool.pkl           # Global compound pool (pickle)
├── pairwise_seqid.tsv          # All-vs-all sequence identity
├── pairwise_pocket_rmsd.tsv    # All-vs-all pocket RMSD
├── passed_targets.txt          # UniProt IDs that passed protein filtering
├── train.txt                   # Train split entries with weights
├── test.txt                    # Test split entries with weights
├── chembl_targets.tsv          # Per-target summary (split, n_actives, n_decoys)
│
└── {UniProt}/                  # Per-target directory
    ├── actives.tsv             # chembl_id, pchembl, smiles
    ├── actives_clustered.tsv   # + cluster_size column
    ├── decoys.tsv              # active_chembl_id -> decoy_ids (;-sep)
    ├── comps/smiles/*.smi      # Raw SMILES files from Stage 1
    ├── pdb/                    # Downloaded PDB + AlphaFold structures
    ├── aligned/                # Structures aligned to AlphaFold model
    ├── pdbid.list              # PDB metadata (method, resolution, chains)
    ├── pocket_info.csv         # Ligand pocket coordinates
    └── sequence.fasta          # Per-target canonical sequence
```

---

## Project Structure

```
ChEMBL-Q/
├── chembl_curator/
│   ├── __init__.py
│   ├── cli.py                  # CLI entry points
│   ├── config.py               # CurationConfig
│   ├── curator.py              # Stage 1
│   ├── protein_filter.py       # Stage 2
│   ├── active_clusterer.py     # Stage 3
│   ├── compound_pool.py        # Stage 4
│   ├── receptor_similarity.py  # Stage 5
│   ├── decoy_selector.py       # Stage 6
│   ├── splitter.py             # Stage 7
│   ├── downloader.py
│   ├── filters.py
│   └── assets/
│       ├── excluded_ligands.txt
│       └── external_targets.fasta
├── docs/
│   └── index.html              # Interactive pipeline overview
├── pyproject.toml
└── README.md
```

---

## External Tools & Databases

- [ChEMBL](https://www.ebi.ac.uk/chembl/) - bioactivity database
- [AlphaFold DB](https://alphafold.ebi.ac.uk/) - predicted protein structures
- [RCSB PDB](https://www.rcsb.org/) - experimental protein structures
- [MMseqs2](https://github.com/soedinglab/MMseqs2) - fast sequence search/clustering
- [nurikit](https://github.com/seoklab/nurikit) - Python TMAlign bindings (used for structure alignment)

## License

This project is provided as-is for research purposes.
