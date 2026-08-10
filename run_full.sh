#!/bin/bash
#SBATCH --job-name=chemblq
#SBATCH -p normal.q
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=96
#SBATCH -o log.out
#SBATCH -e log.err

# Usage:
#   On cluster:  sbatch run_full.sh
#   Locally:     bash run_full.sh
#
# Note: download ChEMBL DB locally before submitting (see CHEMBL_DB below).
#   chembl-curator curate --download --output curated_data_filtered

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────

CONDA_ENV=chemblq
DATA_DIR=curated_data_filtered
N_CPUS=${SLURM_NTASKS:-96}   # auto-read from SLURM; fallback for local runs
N_WORKERS=$N_CPUS
MMseqs_THREADS=$N_CPUS

# Path to pre-downloaded ChEMBL SQLite database
# Download locally first: chembl-curator curate --download --output $DATA_DIR
CHEMBL_DB=/home.galaxy4/j2ho/DB/ChEMBL-Q/chembl_data/chembl_36.db

# External FASTA for co-clustering in stage 7
# Bundled PDBbind+BioLip FASTA is used by default.
# To use a different file: --external-fasta /path/to/your.fasta
# To skip external sequences: --no-external

# ── Environment ───────────────────────────────────────────────────────────────

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "Starting ChEMBL-Q pipeline (${N_CPUS} CPUs, data: ${DATA_DIR})"

# ── Stage 1: Compound filtering  (single-threaded, ~10–30 min) ───────────────

log "Stage 1: compound filtering"
chembl-curator curate \
    --database "$CHEMBL_DB" \
    --config config.json \
    --output "$DATA_DIR" \
    --log-level INFO

# ── Stage 2: Protein filtering  (parallel, ~2–4 h at 96 CPUs) ───────────────

log "Stage 2: protein filtering"
chembl-curator filter-proteins \
    --curated-dir "$DATA_DIR" \
    --n-processes "$N_WORKERS" \
    --log-level INFO

# ── Stage 3: Active clustering  (parallel, ~5 min) ───────────────────────────

log "Stage 3: clustering actives"
chembl-curator cluster-actives \
    --data-dir "$DATA_DIR" \
    --dist-thresh 0.3 \
    --workers "$N_WORKERS" \
    --log-level INFO

# ── Stage 4: Compound pool  (single-threaded, ~1 min) ────────────────────────

log "Stage 4: building compound pool"
chembl-curator build-pool \
    --data-dir "$DATA_DIR" \
    --log-level INFO

# ── Stage 5: Receptor similarity  (parallel, ~30 min) ────────────────────────

log "Stage 5: computing receptor similarity"
chembl-curator receptor-sim \
    --data-dir "$DATA_DIR" \
    --mode both \
    --seqid-threads "$MMseqs_THREADS" \
    --workers "$N_WORKERS" \
    --pocket-radius 10.0 \
    --log-level INFO

# ── Stage 6: Decoy selection  (single-threaded, ~5–15 min) ──────────────────

log "Stage 6: selecting decoys"
chembl-curator select-decoys \
    --data-dir "$DATA_DIR" \
    --max-decoys 30 \
    --seqid-thresh 0.6 \
    --pocket-rmsd-thresh 3.0 \
    --exclusion-mode or \
    --tanimoto-thresh 0.3 \
    --log-level INFO

# ── Stage 7: Train/test split  (MMseqs2, ~5 min) ─────────────────────────────

log "Stage 7: train/test split"
chembl-curator split \
    --data-dir "$DATA_DIR" \
    --seqid 0.3 \
    --valid-frac 0.1 \
    --threads "$MMseqs_THREADS" \
    --log-level INFO

log "Done. Outputs in: $DATA_DIR"
