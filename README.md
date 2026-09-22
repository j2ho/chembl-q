<div align="center">
  <div>&nbsp;</div>
  <img src="docs/chemblq.png" width="300"/>

  <a href="https://j2ho.github.io/chembl-q/">Project Webpage</a>
</div>

A pipeline that curates ChEMBL into a virtual screening benchmark, and a dataset built with it.

Most released scoring models are trained on PDBbind or BioLiP and then scored on DUD-E, DEKOIS 2.0 or LIT-PCBA. Those benchmarks were assembled from targets that have crystal structures, which is the criterion that put the same targets in PDBbind and BioLiP, so the models are being scored on receptors they were trained on. Measured with this pipeline: **99–100% of DUD-E, DEKOIS 2.0 and LIT-PCBA targets have a sequence homologue in PDBbind + BioLiP, and 93–100% have a pocket within 1 Å, at a median distance of 0.00 Å.**

ChEMBL-Q keeps the two apart. 1,081 targets are offered for training and they are the ones those sets already cover; 425 are held back as the benchmark, and none of them has a homologue in PDBbind, in BioLiP, or in the training split itself. Whatever a model was trained on, it can be scored here.

Everything the pipeline measures ships with the data, including the overlap that remains.

---

## The dataset

```
1,506 targets      single-pocket, structure-backed
282,829 actives    exact measurement at or below 10 µM
 60,701 inactives  experimentally measured, kept apart from decoys
   3.1M pairs      property-matched, receptor-aware decoys

train  1,081 targets   97,685 actives   median 17 per target
test     425 targets    8,367 actives   median  2 per target
```

Active counts are clustered representatives, which is what a decoy block is attached to; the unclustered total is 282,829.

**What the test split guarantees.** No test target has a PDBbind or BioLiP homologue at 0.4 identity, none was built from a PDB entry those sets hold, and none is within 0.4 identity of a training sequence. Median distance from a test pocket to its closest external pocket is 2.01 Å, against 0.25 Å for train — a consequence of the sequence rules, since no split decision is made on geometry.

**What it does not guarantee.** The split is on targets, not compounds. 19.1% of test actives are also train actives and 31.1% share a Murcko scaffold with one. Run a ligand-only control before attributing performance to structure; the per-target `actives.tsv` files are there to build a compound-disjoint evaluation from if you need one.

Pocket novelty is reported, not enforced: 188 of the 380 test targets with an external match have a pocket within 2 Å of a PDBbind or BioLiP one, 59 within 1 Å. Filtering on it would remove the data-rich targets, because a fold that has been drugged hard is also one those sets hold many structures of, so `chembl_targets.tsv` carries the per-target distance and leaves the choice to you.

---

## What to score on

The dataset ships as a tarball (**upload it and put the DOI here** — the archive is ~22 MB for `benchmark/` alone, ~49 GB with the structures, so the two are worth publishing separately). `benchmark/` is the part you need to evaluate a model: a labelled screening list per target, with SMILES, so nothing has to be reconstructed and `compound_pool.pkl` never has to be opened.

```
benchmark/targets.tsv           425 rows: counts, pocket overlap, subset flags
benchmark/screen/{UniProt}.tsv  chembl_id, label, smiles, decoy_for

label = active   a binder, exact measurement at or below 10 µM
        decoy    assumed non-binder, an active of a dissimilar receptor
        inactive measured non-binder
```

**Score on the 125 targets with at least five actives** (`pass_reps5`). That is 29% of the held-out targets and 95% of the compounds — 208 of the 425 have exactly one active, where a per-target AUC is either 0 or 1 and averaging them reports noise. All 425 ship anyway: one measured active is still a measurement, and a method that ranks one compound against thirty controls may want them.

|                             | targets | actives | decoys | median actives |
|-----------------------------|--------:|--------:|-------:|---------------:|
| everything held out         |     425 |   8,367 | 247,226 |              2 |
| **at least 5 actives**      | **125** | **7,913** | **234,456** |         **17** |
| + no pocket within 1 Å      |      97 |   6,120 | 180,941 |             17 |
| strict (`pass_strict_2A`)   |      52 |   3,343 |  98,311 |             24 |

The strict set additionally drops every target whose pocket comes within 2 Å of a PDBbind or BioLiP one, and every active that is a training compound or shares its Murcko scaffold with one. Report the primary set; use the strict set to show a result is not an artefact of the overlap that remains.

Two things the table does not say. Decoys and measured inactives are different claims and should not be pooled — a decoy is an assumption, an inactive is a measurement, and there are 7,320 of the latter against 247,226 decoy assignments. And `n_actives_organic_only` counts the actives a conventional featuriser can read: two targets (P0A9P4, P36776) are entirely selenium or boronic-acid chemistry, which MMFF94 has no parameters for. They are recorded rather than removed, because which elements a model can read is the model's property, not the target's.

---

## Pipeline

| Stage | Command | Output |
|-------|---------|--------|
| 1. Compounds and activities | `curate` | `{target}/actives.tsv`, `inactives.tsv`, `measured.tsv` |
| 2. Protein structures | `filter-proteins` | `aligned/*.pdb`, `pocket_info.csv`, `sequences.fasta`, `best_structure.tsv` |
| 3. Active clustering | `cluster-actives` | `{target}/actives_clustered.tsv` |
| 4. Compound pool | `build-pool` | `compound_pool.pkl` |
| 5. Receptor similarity | `receptor-sim` | `pairwise_seqid.tsv`, `pairwise_pocket_hungarian.tsv` |
| 6. Decoy selection | `select-decoys` | `{target}/decoys.tsv` |
| 7. External pocket overlap | `external-pockets`, `pocket-leakage` | `external_pocket_best.tsv`, `external_pocket_hits.tsv` |
| 8. Train/test split | `split` | `train.txt`, `test.txt`, `chembl_targets.tsv` |

Stage 7 writes a table stage 8 reads, so it runs first. Each stage reads the previous stage's output, so a run can be resumed or a single stage re-run with different settings.

---

## Installation

```bash
conda create -n chembl-q python=3.11
conda activate chembl-q
git clone https://github.com/j2ho/chembl-q
cd chembl-q
pip install -e .          # rdkit, numpy, pandas, scipy, requests, tqdm, click, nurikit
pip install -e ".[dev]"   # adds pytest
```

The Python import is `nuri`; the PyPI package is `nurikit`. `pyproject.toml` handles that.

**MMseqs2** must be on `PATH` for stages 5 and 8: [github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2). Structure alignment uses `nurikit`'s TM-align bindings, so no separate TMalign binary is needed. `wget` is used for AlphaFold downloads.

---

## Quick start

```bash
bash run_full.sh          # the whole pipeline at the settings this dataset was built with
```

Or a stage at a time. `config.json` is the scientific contract — without it the defaults are much looser and the result will not match the shipped dataset.

```bash
DATA=curated_v6

chembl-curator curate --database chembl_36.db --config config.json --output $DATA
chembl-curator filter-proteins --curated-dir $DATA --n-processes 32 --cache-dir pdb_cache/
chembl-curator cluster-actives --data-dir $DATA --dist-thresh 0.3 --workers 32
chembl-curator build-pool --data-dir $DATA
chembl-curator receptor-sim --data-dir $DATA --mode seqid --seqid-threads 32
chembl-curator receptor-sim --data-dir $DATA --mode pocket \
    --pocket-method hungarian --pocket-radius 8.0 --workers 32 \
    --output $DATA/pairwise_pocket_hungarian.tsv
chembl-curator select-decoys --data-dir $DATA --max-decoys 30 \
    --pocket-rmsd-thresh 2.0 --min-matched-residues 15 \
    --pocket-rmsd-tsv $DATA/pairwise_pocket_hungarian.tsv
chembl-curator pocket-leakage --data-dir $DATA --cache external_pockets.npz
chembl-curator split --data-dir $DATA --seqid 0.4 --valid-frac 1.0 --threads 32
```

Rough timings at 32 cores: stage 1 ~30 min, stage 2 ~3 h with a warm structure cache, stages 3–5 ~10 min, stage 6 ~45 min, stage 7 ~75 min, stage 8 ~2 min.

---

## Stage reference

### Stage 1: `curate`

Pulls actives, measured inactives and the full tested list out of ChEMBL.

Labelling is by relation, not by threshold alone. An exact measurement (`=`, `<=`) at or below 10 µM is an active. A censored `>` record is non-binding evidence and its strength is the concentration tested, so `>` at or above 100 µM is an inactive. When a compound has both, the conflict is settled by potency rather than by counting records: a compound measured at 0.34 nM across 75 assays stays active despite one depositor calling it inactive, and is tagged contested.

Requiring `pchembl_value >= 5` means only `=` records survive in practice — ChEMBL does not assign a pChEMBL to censored rows.

**Five classes are decided, two are written.** Every compound-target pair is sorted into one of:

| Class | Meaning | Where it lands |
|-------|---------|----------------|
| `active` | exact measurement at or below 10 µM | `actives.tsv` |
| `inactive_annotated` | an explicit inactive `activity_comment` | `inactives.tsv`, evidence `annotated` |
| `inactive_potency` | a censored `>` at or above 100 µM | `inactives.tsv`, evidence `potency` |
| `weak` | exact measurement, but above 10 µM | **nowhere** |
| `uninformative` | only censored records below 100 µM | **nowhere** |

A pair carrying both kinds of inactive evidence is written with evidence `annotated;potency`. In the shipped data: 52,579 `annotated`, 7,651 `potency`, 471 both.

`weak` is the one to know about. Binding was observed and it is real information — the compound is simply too weak to call an active. For a benchmark that is the right thing to drop, because a borderline compound belongs in neither column. For **training** it is a signal being thrown away, and one you cannot recover from the shipped files: a weak compound appears only as a bare ID in `measured.tsv`, indistinguishable from a compound that was tested and found inert. If you want it, re-run stage 1 against `labeler.py`, which already computes the class.

**Verified inactives are not decoys, and the dataset keeps them apart.** Two different kinds of negative ship here:

- `inactives.tsv` — somebody ran the assay and it did not bind. Trustworthy, scarce (60,701 pairs), and biased by what anyone bothered to test.
- `decoys.tsv` — generated: actives of other, dissimilar receptors, property-matched to each active. Plentiful (3.1M pairs) and controllable, but an assumption. Nobody checked that these do not bind.

`measured.tsv` is what keeps the two from mixing. Any compound with *any* ChEMBL record against a target — active, inactive, weak, uninformative — is barred from being that target's decoy, because guessing about a pair someone has direct evidence for is indefensible. The benchmark is scored on the generated decoys; the verified inactives sit alongside so a harder, assumption-free evaluation can be built from them.

```bash
chembl-curator curate --database chembl_36.db --config config.json --output $DATA
chembl-curator curate --create-config myconfig.json   # write an example config
```

Every field of `config.json` is applied; there are no silently-on defaults to remember. The shipped file:

```json
{
  "target_types": ["SINGLE PROTEIN"],
  "activity_types": ["Kd", "Ki", "IC50", "EC50"],
  "relations": ["=", "<="],
  "units": ["nM", "uM"],
  "activity_thresholds": {"nM": 10000.0, "uM": 10.0},
  "min_pchembl_value": 5,
  "min_heavy_atoms": 5,
  "max_heavy_atoms": 80,
  "require_standard_flag": false,
  "exclude_invalid_data": true,
  "exclude_duplicates": true,
  "min_confidence_score": 6,
  "assay_types": ["B"],
  "bao_formats": ["BAO_0000357"],
  "extract_negatives": true,
  "active_max_nm": 10000.0,
  "inactive_min_nm": 100000.0,
  "conflict_decisive_nm": 1000.0
}
```

Actives and negatives share every assay-quality filter. Only `standard_type` differs: an inactive compound has no IC50 to report, so its result is filed under "% Control" or "Inhibition" instead, and restricting negatives to potency types would drop about 90% of depositor inactive calls.

### Stage 2: `filter-proteins`

Fetches structures, aligns them, and keeps targets with a single binding pocket.

```bash
chembl-curator filter-proteins --curated-dir $DATA --n-processes 32 --cache-dir pdb_cache/
```

| Option | Default | Description |
|--------|---------|-------------|
| `--cache-dir` | `pdb_cache/` | Shared structure cache, reused across runs |
| `--use-local-mirror` | off | Try a site-local `pdb_get` before the network |
| `--max-chain-residues` | 1500 | Skip chains longer than this |

RCSB no longer ships PDB format for new entries — 74% of 9-series entries are mmCIF only — so downloads fall back to mmCIF and convert with nuri. The cache matters: this stage deletes the directory of every target it rejects, so without one each run re-downloads several thousand structures it then discards.

A pocket-defining ligand needs at least 5 heavy atoms and a burial of 20 protein atoms within 8 Å, and must not be a modified residue or a membrane glycerophospholipid. Two ligands are in the same pocket when their closest heavy atoms are within 5 Å; centroid distance is size-contaminated and splits elongated cofactors from their own sites. The representative ligand is the one contacting the most residues, ties broken by size, and it is recorded by residue — `Ligand_Residue` in `pocket_info.csv` — not only by its three-letter code. A homo-oligomer holds one copy of the ligand per subunit and the code alone does not say which one the pocket was chosen from.

The receptor chain is bounded from both sides. A UniProt entry frequently maps to a PDB chain holding only a peptide of it, and that peptide is usually the ligand of some other protein in the entry, so `--min-chain-residues` (50) and `--min-chain-coverage` (0.1, as a fraction of the UniProt sequence) reject it. Both are needed: the absolute bar alone passes a 23-residue piece of a 1,426-residue protein, and coverage alone cuts a genuine single domain of a large multi-domain one.

**This stage deletes rejected targets' directories. Never point it at a finished dataset.**

### Stage 3: `cluster-actives`

Butina clustering per target, so a scaffold series counts once. `--dist-thresh 0.3` is Tanimoto similarity ≥ 0.7.

### Stage 4: `build-pool`

Global compound pool: MW, cLogP, TPSA, HBD, HBA, aromatic rings, 2048-bit Morgan fingerprint at radius 2, and the set of targets each compound is active against.

### Stage 5: `receptor-sim`

```bash
chembl-curator receptor-sim --data-dir $DATA --mode seqid --seqid-threads 32
chembl-curator receptor-sim --data-dir $DATA --mode pocket \
    --pocket-method hungarian --pocket-radius 8.0 --workers 32 \
    --output $DATA/pairwise_pocket_hungarian.tsv
```

| Option | Default | Description |
|--------|---------|-------------|
| `--mode` | both | `seqid`, `pocket`, or `both` |
| `--pocket-method` | tmalign | `hungarian` is what the dataset uses |
| `--pocket-radius` | 10.0 | 8.0 with `hungarian` |

A pocket is every residue with a heavy atom within the radius of a ligand heavy atom, reduced to its Cα. `hungarian` pairs residues order-free by assignment on a cost combining a local distance fingerprint with a residue-class penalty, then refines the pairing and the superposition against each other. TM-align's pairing follows sequence order, so two pockets built from the same residues in a different arrangement cannot match: it produced a usable RMSD for 19% of pairs against hungarian's 100%.

### Stage 6: `select-decoys`

```bash
chembl-curator select-decoys --data-dir $DATA --max-decoys 30 \
    --pocket-rmsd-thresh 2.0 --min-matched-residues 15 \
    --pocket-rmsd-tsv $DATA/pairwise_pocket_hungarian.tsv
```

| Option | Default | Description |
|--------|---------|-------------|
| `--max-decoys` | 30 | Decoys per active |
| `--seqid-thresh` | 0.6 | Sequence identity for receptor exclusion |
| `--pocket-rmsd-thresh` | 2.0 | Pocket RMSD for receptor exclusion |
| `--min-matched-residues` | 15 | A close fit over a handful of residues is not evidence |
| `--tanimoto-thresh` | 0.3 | Maximum similarity to the active being paired |
| `--cross-active-thresh` | 0.9 | Maximum similarity to any *other* active of the target |
| `--max-selection-count` | derived | How often one compound may serve as a decoy |

A compound measured against this target is never a decoy for it, whatever the outcome — that is what `measured.tsv` is for. Neither is an active of a target with a similar receptor. Pocket similarity does most of that work: 4,540 pairs qualify on pocket against 240 on sequence.

Property windows: ±50 Da MW, ±2 cLogP, ±50 Å² TPSA, ±2 HBD, ±2 HBA, ±1 aromatic ring.

**Two similarity bars, and they are not equal.** 0.3 applies to the active a decoy is paired with; 0.9 applies to every other active of the same target. Without the second, a decoy for active A was free to be the same molecule as active B, since `measured.tsv` only catches that when the compound happens to have been tested there. The bars cannot both be 0.3: the pool is close to saturated — demand runs at 95.8% of what the reuse cap allows — and refilling a 0.3 cut needs twice the spare capacity, so the run would underfill rather than substitute. 0.9 costs 0.01% of assignments and catches what matters, because the fingerprint is built without chirality and a stereoisomer of a known binder therefore scores 1.0.

`--max-selection-count` is derived as `ceil(total_actives × max_decoys / pool_size) + 1` unless given. Pin it only to reproduce a specific build: a value carried over from a run with a different pool underfills silently. The log reports the value in force, the derived value, and the demand/capacity ratio.

### Stage 7: `external-pockets`, `pocket-leakage`

Superposes every ChEMBL pocket against all 38,825 pockets extracted from PDBbind and BioLiP, using the same pocket definition as stage 5.

```bash
# build the cache once, needs the raw structure sets on disk
chembl-curator external-pockets --biolip-dir /path/BioLiP_updated_set \
    --pdbbind-dir /path/v2020-refined --pdbbind-dir /path/v2020-others \
    --output external_pockets.npz --workers 32

# or use the shipped cache and skip straight to the comparison
chembl-curator pocket-leakage --data-dir $DATA --cache external_pockets.npz --workers 32
```

`external_pockets.npz` is in the repository, so `pocket-leakage` runs without PDBbind or BioLiP on disk.

**The result is reported, not filtered.** Cutting the test set on pocket RMSD removes the data-rich targets: a 1.0 Å cut takes 14% of the test targets but 22% of the actives, because a fold that has been drugged hard is also one those sets hold many structures of. 188 test targets sit within 2 Å of an external pocket and not one of them is the same protein, by construction — the split already guarantees no test target has a homologue at 0.4 identity, so every one of these is a geometric match between unrelated sequences. Stage 8 carries the number into `chembl_targets.tsv` so a benchmark can be scored overall and on a pocket-novel subset.

### Stage 8: `split`

```bash
chembl-curator split --data-dir $DATA --seqid 0.4 --valid-frac 1.0 --threads 32
chembl-curator split --data-dir $DATA --no-external          # ChEMBL-only split
chembl-curator split --data-dir $DATA --external-fasta my.fa # your own reference set
```

A target is kept out of the test set if any of three rules fires. All three are sequence or identifier rules; pocket geometry is never used to decide a split.

1. **Sequence homology**, 978 targets. An external entry aligns at 0.4 identity or better over at least 80% of that external sequence. Coverage is measured on the external side because what gets crystallised is a domain while a ChEMBL target is a whole UniProt entry; measuring over the query let 120 targets into the test set with their own PDB entry sitting in BioLiP.
2. **Shared structure**, 680 targets. The target was built from a PDB entry those sets hold. Alignment cannot be relied on here: MMseqs2 finds only a 31-residue alignment between IGF2R's 2,491 residues and its 182-residue PDBbind construct.
3. **Proximity to train**, 3 targets. Iterated to a fixed point.

`--valid-frac 1.0` sends every eligible cluster to test: the test set is defined by what is safe to hold out, not by a target ratio.

**External FASTA IDs** are dot-prefixed as `>{source}.{entry_id}`:

```
>pdbbind.1a4k
>biolip.10gs_VWW_A_1
>myscreendb.custom_entry_42
```

---

## Output

```
curated_v6/
├── sequences.fasta                # canonical sequences, all passed targets
├── best_structure.tsv             # uniprot -> best ligand-bound chain + resolution
├── passed_targets.txt
├── compound_pool.pkl
├── pairwise_seqid.tsv             # MMseqs2 all-vs-all
├── pairwise_pocket_hungarian.tsv  # 1,155,960 internal pocket pairs
├── external_pocket_best.tsv       # closest PDBbind/BioLiP pocket per target
├── external_pocket_hits.tsv       # every pair under the reported cutoff
├── train.txt / test.txt           # one line per active, with a sampling weight
├── chembl_targets.tsv             # split, counts, and the external pocket match
├── benchmark/                     # the released evaluation set, see above
│   ├── targets.tsv                #   425 rows, subset flags
│   ├── screen/{UniProt}.tsv       #   chembl_id, label, smiles, decoy_for
│   └── README.md
│
└── {UniProt}/
    ├── actives.tsv                # chembl_id, pchembl, smiles
    ├── actives_clustered.tsv      # + cluster_size
    ├── inactives.tsv              # chembl_id, evidence, smiles — measured, not generated
    │                              #   evidence: annotated | potency | annotated;potency
    ├── measured.tsv               # every compound tested here, whatever the label;
    │                              #   none of these can be this target's decoy
    ├── decoys.tsv                 # active_chembl_id -> decoy_ids (;-separated)
    ├── pocket_info.csv            # structure, chain, ligand code AND residue, centre
    ├── aligned/                   # superposed chain PDBs
    ├── pdb/                       # downloaded structures
    ├── pdbid.list
    └── sequence.fasta
```

`chembl_targets.tsv` gains two columns when stage 7 has run, and is omitted entirely when it has not — a 0.0 default would read as "identical to an external pocket":

```
uniprot  split  n_actives  n_decoys  ext_pocket_rmsd  ext_pocket_entry
Q9Y5N1   test         409     12265            2.534  biolip.4qkx_35V_A_1
```

---

## Project structure

```
ChEMBL-Q/
├── chembl_curator/
│   ├── cli.py                  # CLI entry points
│   ├── config.py               # CurationConfig
│   ├── curator.py              # stage 1
│   ├── labeler.py              # active / inactive / contested, no DB access
│   ├── protein_filter.py       # stage 2
│   ├── active_clusterer.py     # stage 3
│   ├── compound_pool.py        # stage 4
│   ├── receptor_similarity.py  # stage 5
│   ├── pocket_align.py         # order-free pocket superposition
│   ├── decoy_selector.py       # stage 6
│   ├── external_pockets.py     # stage 7
│   ├── splitter.py             # stage 8
│   └── assets/
│       ├── excluded_ligands.txt
│       └── external_targets.fasta
├── tests/                      # pytest, no network needed except the download tests
├── config.json
├── external_pockets.npz        # 38,825 PDBbind/BioLiP pockets, 17 MB
├── run_full.sh
└── docs/index.html
```

```bash
pytest tests/ -q
```

---

## Built on

- [ChEMBL](https://www.ebi.ac.uk/chembl/) — bioactivity data
- [RCSB PDB](https://www.rcsb.org/) and [AlphaFold DB](https://alphafold.ebi.ac.uk/) — structures
- [MMseqs2](https://github.com/soedinglab/MMseqs2) — sequence search and clustering
- [nurikit](https://github.com/seoklab/nurikit) — molecule I/O and TM-align bindings
- [RDKit](https://www.rdkit.org/) — cheminformatics
- PDBbind v2020 and BioLiP — the reference sets held out against

## License

Provided as-is for research purposes.
