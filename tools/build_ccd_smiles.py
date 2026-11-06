#!/usr/bin/env python
"""Build chembl_curator/assets/ccd_smiles.tsv.gz from the PDB chemical
component dictionary.

A cocrystal ligand is recorded in pocket_info.csv only as a three-letter PDB
component id, which carries no structure. Stage 4 needs SMILES for it so the
cross-active bar can see the one binder of a target that is certain: the
molecule a crystallographer found in the pocket we define.

The CCD is ~120 MB compressed and is NOT vendored. This script distills it to
the two columns the pipeline uses, which gzip to about 1 MB.

    python tools/build_ccd_smiles.py                    # download, then write
    python tools/build_ccd_smiles.py components.cif.gz  # from a local copy

Re-run it to pick up components added to the PDB since the asset was built.
"""

import argparse
import gzip
import re
import sys
import urllib.request
from pathlib import Path

CCD_URL = "https://files.rcsb.org/pub/pdb/data/monomers/components.cif.gz"
ASSET = (Path(__file__).resolve().parent.parent / "chembl_curator" /
         "assets" / "ccd_smiles.tsv.gz")

# comp_id  type  program  version  descriptor
ROW = re.compile(r'^(\S+)\s+(SMILES_CANONICAL|SMILES)\s+(\S+)\s+(\S+)\s+(.*)$')
PREFERRED_PROGRAMS = ('OpenEye', 'CACTVS')

# The CCD held 52,070 components with a usable SMILES when this was written.
# A parse that collapses far below that means the format moved, and writing the
# asset anyway would quietly shrink the cross-active reference set.
MIN_COMPONENTS = 30_000


def parse(src: Path) -> dict:
    """comp_id -> SMILES, preferring the canonical OpenEye/CACTVS form.

    The CCD carries several SMILES per component from different toolkits, and
    they are not interchangeable: the ACDLabs rows often drop stereochemistry.
    Rank the candidates rather than taking whichever comes first.
    """
    best: dict = {}
    with gzip.open(src, 'rt', errors='replace') as handle:
        for line in handle:
            if ' SMILES' not in line:
                continue
            match = ROW.match(line.strip())
            if not match:
                continue
            comp_id, kind, program, _version, desc = match.groups()
            desc = desc.strip()
            if desc[:1] in ('"', "'") and desc[-1:] == desc[:1]:
                desc = desc[1:-1]
            if not desc or desc == '?' or '\t' in desc:
                continue
            rank = (0 if kind == 'SMILES_CANONICAL' else 1,
                    0 if program in PREFERRED_PROGRAMS else 1)
            previous = best.get(comp_id)
            if previous is None or rank < previous[0]:
                best[comp_id] = (rank, desc)
    return {k: v[1] for k, v in best.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', nargs='?', type=Path,
                        help='components.cif.gz; downloaded if omitted')
    parser.add_argument('--output', type=Path, default=ASSET)
    args = parser.parse_args()

    source = args.source
    downloaded = None
    if source is None:
        downloaded = args.output.parent / 'components.cif.gz.tmp'
        downloaded.parent.mkdir(parents=True, exist_ok=True)
        print(f'downloading {CCD_URL}')
        urllib.request.urlretrieve(CCD_URL, downloaded)
        source = downloaded

    smiles = parse(source)
    if len(smiles) < MIN_COMPONENTS:
        print(f'only {len(smiles)} components parsed, refusing to write the '
              f'asset: expected at least {MIN_COMPONENTS}', file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, 'wt') as handle:
        handle.write('comp_id\tsmiles\n')
        for comp_id in sorted(smiles):
            handle.write(f'{comp_id}\t{smiles[comp_id]}\n')
    if downloaded is not None:
        downloaded.unlink()

    size_mb = args.output.stat().st_size / 1e6
    print(f'{len(smiles)} components -> {args.output} ({size_mb:.1f} MB)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
