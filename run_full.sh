#!/bin/bash
#SBATCH --job-name=chemblq
#SBATCH -p normal.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH -o log.out
#SBATCH -e log.err

# End-to-end ChEMBL-Q build, at the settings the released dataset was built
# with. Every parameter here is deliberate; see README.md for what each one
# does and why.
#
# Usage:
#   sbatch run_full.sh          # cluster
#   bash run_full.sh            # local
#
# Requirements beyond the Python package:
#   mmseqs on PATH                      (stages 5 and 8)
#   the ChEMBL SQLite database          (stage 1, see CHEMBL_DB)
#   PDBbind and BioLiP raw structures   (stage 7a only; skip it and stage 7b
#                                        if you do not have them, the split
#                                        then simply omits two columns)

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────

CONDA_ENV=chembl-q
DATA_DIR=${DATA_DIR:-curated_v5}
N_CPUS=${SLURM_CPUS_PER_TASK:-8}

# Pre-downloaded ChEMBL SQLite. To fetch it:
#   chembl-curator curate --download --output "$DATA_DIR"
CHEMBL_DB=${CHEMBL_DB:-chembl_data/chembl_36.db}

# Structure cache, shared across runs. Stage 2 discards the directory of every
# target it rejects, so without this each run re-downloads ~4,400 structures it
# then throws away; most of a 5 h stage 2 was network.
CACHE_DIR=${CACHE_DIR:-pdb_cache}

# Raw external structure sets, for stage 7a. Leave empty to skip stages 7a/7b.
BIOLIP_DIR=${BIOLIP_DIR:-}
PDBBIND_REFINED=${PDBBIND_REFINED:-}
PDBBIND_OTHERS=${PDBBIND_OTHERS:-}
EXT_POCKETS=${EXT_POCKETS:-external_pockets.npz}

# ── Environment ───────────────────────────────────────────────────────────────

if command -v conda >/dev/null; then
    . "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
log "ChEMBL-Q build: ${N_CPUS} CPUs, data dir ${DATA_DIR}"

# ── Stage 1: compound and activity curation  (~30 min) ───────────────────────
# Actives, experimental inactives, and the measured/conflict tables. The config
# is the scientific contract: BAO_0000357, assay type B, confidence >= 6,
# pchembl >= 5. Running without --config uses much looser defaults and will not
# reproduce the released dataset.

log "Stage 1: curating compounds and activities"
chembl-curator curate \
    --database "$CHEMBL_DB" \
    --config config.json \
    --output "$DATA_DIR" \
    --log-level INFO

# ── Stage 2: protein structures  (~3 h at 16 CPUs, cache warm) ──────────────
# Downloads structures (mmCIF fallback included), aligns them, keeps targets
# with a single binding pocket, and picks one representative ligand per target.
# THIS DELETES THE DIRECTORY OF EVERY REJECTED TARGET, so never point it at a
# finished dataset.

log "Stage 2: protein structure filtering"
chembl-curator filter-proteins \
    --curated-dir "$DATA_DIR" \
    --n-processes "$N_CPUS" \
    --cache-dir "$CACHE_DIR" \
    --log-level INFO

# ── Stage 3: active clustering  (~1 min) ─────────────────────────────────────

log "Stage 3: clustering actives"
chembl-curator cluster-actives \
    --data-dir "$DATA_DIR" \
    --dist-thresh 0.3 \
    --workers "$N_CPUS" \
    --log-level INFO

# ── Stage 4: compound pool  (~1 min) ─────────────────────────────────────────

log "Stage 4: building compound pool"
chembl-curator build-pool \
    --data-dir "$DATA_DIR" \
    --log-level INFO

# ── Stage 5: receptor similarity  (~6 min) ───────────────────────────────────
# Sequence identity, then order-free pocket RMSD. The Hungarian method at an
# 8 A heavy-atom pocket radius is what the dataset uses: TM-align produced a
# usable RMSD for only 19% of pairs because its residue pairing follows
# sequence order, and cross-fold pocket similarity is the case that matters.

log "Stage 5a: sequence identity"
chembl-curator receptor-sim \
    --data-dir "$DATA_DIR" --mode seqid \
    --seqid-threads "$N_CPUS" --log-level INFO

log "Stage 5b: pocket RMSD"
chembl-curator receptor-sim \
    --data-dir "$DATA_DIR" --mode pocket \
    --pocket-method hungarian --pocket-radius 8.0 \
    --workers "$N_CPUS" \
    --output "$DATA_DIR/pairwise_pocket_hungarian.tsv" \
    --log-level INFO

# ── Stage 6: decoy selection  (~45 min) ──────────────────────────────────────
# A compound measured against this target is never a decoy for it, and neither
# is an active of a target with a similar receptor. --min-matched-residues 15
# rejects a close fit that rests on a handful of residues.

log "Stage 6: selecting decoys"
chembl-curator select-decoys \
    --data-dir "$DATA_DIR" \
    --max-decoys 30 \
    --seqid-thresh 0.6 \
    --pocket-rmsd-thresh 2.0 \
    --min-matched-residues 15 \
    --pocket-rmsd-tsv "$DATA_DIR/pairwise_pocket_hungarian.tsv" \
    --exclusion-mode or \
    --tanimoto-thresh 0.3 \
    --seed 42 \
    --log-level INFO

# ── Stage 7: pocket-level overlap with PDBbind and BioLiP  (~75 min) ────────
# Measured and reported, not filtered on: cutting the test set by pocket RMSD
# removes the data-rich targets, because a fold that has been drugged hard is
# also one those sets hold many structures of. Stage 8 carries the per-target
# result into chembl_targets.tsv so a benchmark can be scored overall and on a
# pocket-novel subset.

if [ -n "$BIOLIP_DIR" ] && [ -n "$PDBBIND_REFINED" ]; then
    if [ ! -f "$EXT_POCKETS" ]; then
        log "Stage 7a: building the external pocket cache"
        chembl-curator external-pockets \
            --biolip-dir "$BIOLIP_DIR" \
            --pdbbind-dir "$PDBBIND_REFINED" \
            --pdbbind-dir "$PDBBIND_OTHERS" \
            --output "$EXT_POCKETS" \
            --pocket-radius 8.0 --workers "$N_CPUS" --log-level INFO
    else
        log "Stage 7a: reusing $EXT_POCKETS"
    fi

    log "Stage 7b: scoring ChEMBL pockets against the external set"
    chembl-curator pocket-leakage \
        --data-dir "$DATA_DIR" \
        --cache "$EXT_POCKETS" \
        --pocket-radius 8.0 \
        --rmsd-report 2.0 \
        --min-matched-residues 15 \
        --workers "$N_CPUS" --log-level INFO
else
    log "Stage 7: skipped, no BIOLIP_DIR/PDBBIND_REFINED set"
fi

# ── Stage 8: train/test split  (~2 min) ──────────────────────────────────────
# Test is ChEMBL-only and separated from PDBbind and BioLiP by sequence.
# valid-frac 1.0 sends every eligible cluster to test; the test set is defined
# by what is safe to hold out, not by a target ratio.
# Reads stage 7's output when present, and must therefore run after it.

log "Stage 8: train/test split"
chembl-curator split \
    --data-dir "$DATA_DIR" \
    --seqid 0.4 \
    --valid-frac 1.0 \
    --threads "$N_CPUS" \
    --log-level INFO

log "Done. Outputs in: $DATA_DIR"
