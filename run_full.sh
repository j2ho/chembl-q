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
# normal.q mixes Intel Ice Lake (AVX512) and AMD Zen 2 (no AVX512), and the
# local mmseqs build needs AVX512: on an AMD node it dies with Illegal
# instruction, taking stage 5a and stage 8 with it. SLURM defines no node
# features here, so --exclude is the only lever. star039-046 and star050 were
# never sampled and are excluded defensively.
#SBATCH --exclude=star019,star021,star023,star025,star027,star030,star032,star034,star035,star037,star038,star039,star040,star041,star042,star043,star044,star045,star047,star048,star049,star050

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
DATA_DIR=${DATA_DIR:-chembl-q-data}
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

# Python salts string hashing per process, so iterating any set of compound or
# target IDs comes out in a different order every run. Stage 6 draws decoys
# against a running usage cap, which makes that order change which decoys get
# picked; seeding random covers the shuffle but not the iteration. The code
# sorts where order decides an output, and this pins the rest.
export PYTHONHASHSEED=0

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

# ── Stage 4: compound pool  (~2 min) ─────────────────────────────────────────
# Also builds the cross-active reference set stage 6 compares against, which
# needs weak.tsv from stage 1 and pocket_info.csv from stage 2. Re-run it after
# either of those, or stage 6 refuses to start.

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
#
# --max-selection-count is left derived on purpose. It is tight enough to
# matter: on the released build the derived cap was 42, demand ran at 96% of
# capacity, and most of the pool sat exactly on the cap, so a value pinned
# from a run with a different pool underfills silently rather than erroring.
# The run logs the value in force, the derived value, and the demand/capacity
# ratio. Pin it here only to reproduce one specific build.
#
# --cross-active-thresh 0.5 is the bar against the target's whole reference
# set, separate from the 0.3 against the paired active. Stage 4 builds that
# reference from every active (not the cluster representatives), every
# same-pocket cocrystal ligand, and every weak binder; measured inactives are
# deliberately left out, since resembling a confirmed non-binder is a point in
# a decoy's favour.
#
# Measured on the previous build, against the full reference: 16.34% of
# assignments sat above 0.3, 0.855% above 0.5, 0.097% above 0.7. 0.3 needs
# twice the spare capacity the reuse cap leaves and ends in underfill; 0.5
# refilled to 3,140,664 pairs with only 2,016 actives short.
#
# The reference set matters more than the number. At 0.9 against
# representatives only, 94 assignments still reached above 0.9 similarity to a
# non-representative active of their own target - one of them an exact
# duplicate - and 14 were the same molecule as a cocrystal ligand of the very
# target they were a decoy for.

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
    --cross-active-thresh 0.5 \
    --seed 42 \
    --log-level INFO

# ── Stage 7: pocket-level overlap with PDBbind and BioLiP  (~75 min) ────────
# Measured and reported, not filtered on: cutting the test set by pocket RMSD
# removes the data-rich targets, because a fold that has been drugged hard is
# also one those sets hold many structures of. Stage 8 carries the per-target
# result into chembl_targets.tsv so a benchmark can be scored overall and on a
# pocket-novel subset.

# The cache is built from external structures alone and carries no ChEMBL
# content, so it survives a rebuild of this dataset. Only 7a needs the raw
# directories; 7b needs the cache. Requiring the directories for both used to
# skip the whole stage on a rerun that had the cache but not the sources.
if [ ! -f "$EXT_POCKETS" ]; then
    if [ -n "$BIOLIP_DIR" ] && [ -n "$PDBBIND_REFINED" ]; then
        log "Stage 7a: building the external pocket cache"
        chembl-curator external-pockets \
            --biolip-dir "$BIOLIP_DIR" \
            --pdbbind-dir "$PDBBIND_REFINED" \
            --pdbbind-dir "$PDBBIND_OTHERS" \
            --output "$EXT_POCKETS" \
            --pocket-radius 8.0 --workers "$N_CPUS" --log-level INFO
    else
        log "Stage 7: skipped, no $EXT_POCKETS and no BIOLIP_DIR/PDBBIND_REFINED"
    fi
else
    log "Stage 7a: reusing $EXT_POCKETS"
fi

if [ -f "$EXT_POCKETS" ]; then
    log "Stage 7b: scoring ChEMBL pockets against the external set"
    chembl-curator pocket-leakage \
        --data-dir "$DATA_DIR" \
        --cache "$EXT_POCKETS" \
        --pocket-radius 8.0 \
        --rmsd-report 2.0 \
        --min-matched-residues 15 \
        --workers "$N_CPUS" --log-level INFO
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
