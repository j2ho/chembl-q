#!/usr/bin/env python
"""Write a release's benchmark/ directory: targets.tsv and screen/*.tsv.

The script that produced v6's benchmark was never committed, so this is a
reconstruction, and it is validated the only way a reconstruction can be: by
regenerating v6 and diffing against the shipped file.

    python tools/build_benchmark.py --data curated_v6 --out /tmp/bench_v6
    diff curated_v6/benchmark/targets.tsv /tmp/bench_v6/targets.tsv

Two definitions were recovered that way rather than guessed, and both matter:

  * n_actives_no_train_compound counts against EVERY active of a training
    target, not the cluster representatives. On Q9Y5N1 the representative set
    gives 354 and the full set gives the shipped 349.
  * n_decoys is counted from decoys.tsv, not derived as n_actives * 30. v6's
    test split holds 247,226 assignments where the derived figure would claim
    251,010; the gap is actives the pool underfilled.
  * pass_strict_2A is NOT pass_reps5 AND pass_pocket_2A, which is the obvious
    reading and gives 66 targets. It asks for 5 actives whose Murcko scaffold
    is absent from the training split, which gives the shipped 52. Requiring
    5 novel COMPOUNDS as well would also reproduce 52, but novel-scaffold
    implies novel-compound (a compound present in training carries its own
    scaffold there, and the data confirms it: no target has more novel
    scaffolds than novel compounds), so that conjunct is dead weight.
  * a ringless active has an EMPTY Murcko scaffold, which is a member of
    every training set, so its novelty falls back to the compound id.
    See is_scaffold_novel.

Only test-split targets ship. `train.txt` and `test.txt` are per-active lists
over three sources, so the target list is the distinct entry_id of their
`chembl` rows.
"""

import argparse
import csv
import pickle
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog('rdApp.*')

# What n_actives_organic_only counts. B and Se are what this excludes in
# practice; they are also what a MotifScreen featuriser has no atom type for.
ORGANIC = {'H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'}

COLUMNS = (
    'uniprot', 'n_actives', 'n_decoys', 'n_decoys_unique', 'n_inactives',
    'ext_pocket_rmsd', 'ext_pocket_entry', 'n_actives_organic_only',
    'n_actives_no_train_compound', 'n_actives_no_train_scaffold',
    'pass_reps5', 'pass_pocket_2A', 'pass_pocket_1A', 'pass_strict_2A',
)

# The bar for both pass_reps5 (actives) and pass_strict_2A (novel-scaffold
# actives). A per-target AUC over fewer than five positives is noise.
MIN_REPS = 5
SCREEN_HEADER = ('chembl_id', 'label', 'smiles', 'decoy_for')


def split_targets(path: Path, source: str = 'chembl') -> list:
    """Distinct entry_id of one source's rows, sorted."""
    out = set()
    with open(path) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0] == source:
                out.add(parts[1])
    return sorted(out)


def read_pairs(path: Path, smiles_column: str = 'smiles') -> list:
    """[(chembl_id, smiles)] of a headed TSV, file order, duplicates kept."""
    rows = []
    if not path.exists():
        return rows
    with open(path) as handle:
        header = next(handle, '').rstrip('\n').split('\t')
        if smiles_column not in header:
            return rows
        at = header.index(smiles_column)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) > at:
                rows.append((parts[0], parts[at]))
    return rows


def scaffold_of(smiles: str):
    """Murcko scaffold SMILES, or None when the molecule has no ring system.

    None is not "failed to compute": a peptide or a polyamine has an EMPTY
    Murcko scaffold, and an empty scaffold is in every training set, so a
    membership test on it is vacuously true. Callers fall back to compound
    identity instead.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    core = MurckoScaffold.GetScaffoldForMol(mol)
    if core is None:
        return None
    return Chem.MolToSmiles(core) or None


def is_scaffold_novel(cid: str, smiles: str, compounds: set,
                      scaffolds: set) -> bool:
    """Scaffold unseen in training, falling back to the compound when ringless.

    P42575 is the one target in v6's 425 where this matters: three of its 21
    actives are ringless (two caspase substrate peptides and spermine). Two are
    training compounds and one is not, which is how the shipped 13 sits between
    the 12 you get calling them all seen and the 15 you get calling them all
    novel.
    """
    scaffold = scaffold_of(smiles)
    if scaffold is None:
        return cid not in compounds
    return scaffold not in scaffolds


def is_organic(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return all(atom.GetSymbol() in ORGANIC for atom in mol.GetAtoms())


def train_reference(data: Path, train: list):
    """Every active compound id and Murcko scaffold of the training split."""
    compounds, scaffolds = set(), set()
    for i, uniprot in enumerate(train, 1):
        for cid, smiles in read_pairs(data / uniprot / 'actives.tsv'):
            if cid in compounds:
                continue
            compounds.add(cid)
            scaffold = scaffold_of(smiles)
            if scaffold:
                scaffolds.add(scaffold)
        if i % 100 == 0:
            print(f'  train {i}/{len(train)}: {len(compounds)} compounds, '
                  f'{len(scaffolds)} scaffolds', flush=True)
    return compounds, scaffolds


def external_pockets(path: Path) -> dict:
    out = {}
    if not path.exists():
        return out
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            out[row['chembl_target']] = (row['pocket_rmsd'],
                                         row['external_entry'])
    return out


def decoys_of(path: Path) -> dict:
    """active -> [decoy ids], in the order stage 6 wrote them."""
    out = {}
    if not path.exists():
        return out
    with open(path) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                out[parts[0]] = [d for d in parts[1].split(';') if d]
    return out


def pool_smiles(path: Path) -> dict:
    """chembl_id -> SMILES for every poolable compound.

    decoys.tsv stores ids only, and a decoy is by construction a cluster
    representative of some target, so the stage 4 pool is the one place that
    certainly has its structure. Reading 1,506 actives files instead would
    work but would reconstruct what the pickle already holds.
    """
    with open(path, 'rb') as handle:
        pool = pickle.load(handle)['pool']
    return {cid: rec['smiles'] for cid, rec in pool.items()}


def write_screen(out_dir: Path, uniprot: str, actives: list, decoys: dict,
                 inactives: list, smiles: dict) -> Counter:
    """One row per labelled compound; a decoy repeats once per active it serves.

    The repetition is the point: `decoy_for` makes the per-active comparison
    reachable, and 30 controlled comparisons per active is the unit the
    benchmark reports. A consumer that wants the target-level list drops the
    column and deduplicates.
    """
    counts = Counter()
    active_smiles = dict(actives)
    path = out_dir / f'{uniprot}.tsv'
    with open(path, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(SCREEN_HEADER)
        for cid, smi in actives:
            writer.writerow([cid, 'active', smi, ''])
            counts['active'] += 1
        for active_id in sorted(decoys):
            for cid in decoys[active_id]:
                smi = smiles.get(cid) or active_smiles.get(cid)
                if not smi:
                    counts['decoy_no_smiles'] += 1
                    continue
                writer.writerow([cid, 'decoy', smi, active_id])
                counts['decoy'] += 1
        for cid, smi in inactives:
            writer.writerow([cid, 'inactive', smi, ''])
            counts['inactive'] += 1
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--targets-only', action='store_true',
                        help='write targets.tsv and skip screen/')
    args = parser.parse_args()

    data = args.data
    test = split_targets(data / 'test.txt')
    train = split_targets(data / 'train.txt')
    if not test or not train:
        print(f'empty split under {data}', file=sys.stderr)
        return 1
    print(f'{len(test)} test targets, {len(train)} train targets', flush=True)

    ext = external_pockets(data / 'external_pocket_best.tsv')
    compounds, scaffolds = train_reference(data, train)
    print(f'train reference: {len(compounds)} compounds, {len(scaffolds)} '
          f'scaffolds', flush=True)

    smiles = {} if args.targets_only else pool_smiles(
        data / 'compound_pool.pkl')
    args.out.mkdir(parents=True, exist_ok=True)
    screen_dir = args.out / 'screen'
    if not args.targets_only:
        screen_dir.mkdir(exist_ok=True)

    rows, totals = [], Counter()
    for i, uniprot in enumerate(test, 1):
        tdir = data / uniprot
        actives = read_pairs(tdir / 'actives_clustered.tsv')
        decoys = decoys_of(tdir / 'decoys.tsv')
        inactives = read_pairs(tdir / 'inactives.tsv')

        unique = set()
        for assigned in decoys.values():
            unique.update(assigned)

        no_scaffold = sum(1 for cid, smi in actives
                          if is_scaffold_novel(cid, smi, compounds, scaffolds))

        rmsd, entry = ext.get(uniprot, ('', ''))
        value = float(rmsd) if rmsd else None
        reps5 = len(actives) >= MIN_REPS
        pocket2 = value is None or value > 2.0
        pocket1 = value is None or value > 1.0
        strict = no_scaffold >= MIN_REPS and pocket2
        rows.append({
            'uniprot': uniprot,
            'n_actives': len(actives),
            'n_decoys': sum(len(v) for v in decoys.values()),
            'n_decoys_unique': len(unique),
            'n_inactives': len(inactives),
            'ext_pocket_rmsd': rmsd,
            'ext_pocket_entry': entry,
            'n_actives_organic_only': sum(1 for _, s in actives
                                          if is_organic(s)),
            'n_actives_no_train_compound': sum(1 for cid, _ in actives
                                               if cid not in compounds),
            'n_actives_no_train_scaffold': no_scaffold,
            'pass_reps5': int(reps5),
            'pass_pocket_2A': int(pocket2),
            'pass_pocket_1A': int(pocket1),
            'pass_strict_2A': int(strict),
        })
        if not args.targets_only:
            totals += write_screen(screen_dir, uniprot, actives, decoys,
                                   inactives, smiles)
        if i % 50 == 0:
            print(f'  target {i}/{len(test)}', flush=True)

    with open(args.out / 'targets.tsv', 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(COLUMNS),
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

    flags = {name: sum(r[name] for r in rows) for name in
             ('pass_reps5', 'pass_pocket_2A', 'pass_pocket_1A',
              'pass_strict_2A')}
    print(f'wrote {len(rows)} rows to {args.out / "targets.tsv"}')
    print(f'  flags: {flags}')
    if not args.targets_only:
        print(f'  screen rows: {dict(totals)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
