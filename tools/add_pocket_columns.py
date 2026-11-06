#!/usr/bin/env python
"""Add the pocket-flag columns to chembl_targets.tsv, for all 1,506 targets.

benchmark/targets.tsv covers the 425 test targets only, and the targets whose
pocket sits on a cofactor are mostly in TRAINING: 3,294 train rows against 81
test rows for the displaced class, 806 against 124 for the unverifiable one.
So a training loader had no way to see the flag. This writes it where every
target is listed.

Derived entirely from pocket_info.csv, so no pocket, decoy or split is
recomputed and every other column is copied through unchanged.

    python tools/add_pocket_columns.py --data curated_v7
    python tools/add_pocket_columns.py --data curated_v7 --check
"""

import argparse
import shutil
import sys
from pathlib import Path

from chembl_curator.pocket_flags import (
    pass_pocket_noncofactor,
    pocket_cofactor_status,
    pocket_ligands,
    prosthetic_codes,
)

NEW_COLUMNS = ('pocket_ligand', 'pass_pocket_noncofactor',
               'pocket_cofactor_status')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--check', action='store_true',
                        help='report the counts and write nothing')
    args = parser.parse_args()

    path = args.data / 'chembl_targets.tsv'
    if not path.exists():
        print(f'no such file: {path}', file=sys.stderr)
        return 1
    prosthetic = prosthetic_codes()
    if not prosthetic:
        print('prosthetic group asset missing; refusing to write blank '
              'columns', file=sys.stderr)
        return 1

    with open(path) as handle:
        lines = [l.rstrip('\n') for l in handle if l.strip()]
    header = lines[0].split('\t')
    already = [c for c in NEW_COLUMNS if c in header]
    if already and not args.check:
        print(f'{path} already carries {already}; nothing to do')
        return 0

    out = ['\t'.join(header + list(NEW_COLUMNS))]
    counts = {'displaced': 0, 'unverifiable': 0, '': 0}
    by_split = {}
    for line in lines[1:]:
        parts = line.split('\t')
        uniprot = parts[0]
        split = parts[1] if len(parts) > 1 else '?'
        tdir = args.data / uniprot
        status = pocket_cofactor_status(tdir, prosthetic)
        counts[status] = counts.get(status, 0) + 1
        if status:
            by_split.setdefault(status, {}).setdefault(split, 0)
            by_split[status][split] += 1
        out.append('\t'.join(parts + [
            ';'.join(pocket_ligands(tdir)),
            str(pass_pocket_noncofactor(status)),
            status,
        ]))

    print(f'{len(lines) - 1} targets')
    for status in ('displaced', 'unverifiable'):
        print(f'  {status:13s} {counts.get(status, 0):4d}   '
              f'{by_split.get(status, {})}')
    print(f'  ordinary      {counts.get("", 0):4d}')

    if args.check:
        print('--check: nothing written')
        return 0
    backup = path.with_suffix('.tsv.before_pocket_columns')
    if not backup.exists():
        shutil.copy2(path, backup)
    path.write_text('\n'.join(out) + '\n')
    print(f'wrote {path} ({len(NEW_COLUMNS)} columns added; original kept at '
          f'{backup.name})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
