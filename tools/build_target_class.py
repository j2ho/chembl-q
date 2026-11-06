#!/usr/bin/env python
"""Build chembl_curator/assets/target_class.tsv.gz from the ChEMBL database.

The benchmark reports n_decoys_same_family: assignments where the decoy is a
documented active of a target sharing this one's protein classification path.
That needs a uniprot -> classification map, which lives in a 30 GB SQLite file
nobody downstream should have to hold. The map itself is a few hundred KB.

    python tools/build_target_class.py chembl_data/chembl_36.db

One row per (accession, path) pair, since a protein can sit under several
classification paths and the comparison is set against set. Paths contain
spaces but never tabs, so no escaping is needed.

Re-run it against a newer ChEMBL to refresh the asset.
"""

import argparse
import gzip
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ASSET = (Path(__file__).resolve().parent.parent / "chembl_curator" /
         "assets" / "target_class.tsv.gz")

QUERY = """
SELECT cs.accession, pc.protein_class_desc
FROM component_sequences cs
JOIN component_class cc ON cs.component_id = cc.component_id
JOIN protein_classification pc ON cc.protein_class_id = pc.protein_class_id
WHERE cs.accession IS NOT NULL AND pc.protein_class_desc IS NOT NULL
"""

# ChEMBL 36 classified 12,818 accessions when this was written. A collapse far
# below that means the schema moved, and writing the asset anyway would report
# a family overlap of zero and look like good news.
MIN_ACCESSIONS = 8_000


def read_db(path: Path) -> dict:
    out = defaultdict(set)
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        for accession, desc in con.execute(QUERY):
            desc = desc.strip()
            if desc and "\t" not in desc:
                out[accession].add(desc)
    finally:
        con.close()
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("--output", type=Path, default=ASSET)
    args = parser.parse_args()

    if not args.database.exists():
        print(f"no such database: {args.database}", file=sys.stderr)
        return 1
    classes = read_db(args.database)
    if len(classes) < MIN_ACCESSIONS:
        print(f"only {len(classes)} accessions classified, refusing to write "
              f"the asset: expected at least {MIN_ACCESSIONS}", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    with gzip.open(args.output, "wt") as handle:
        handle.write("accession\tprotein_class_desc\n")
        for accession in sorted(classes):
            for desc in sorted(classes[accession]):
                handle.write(f"{accession}\t{desc}\n")
                n_rows += 1

    size_kb = args.output.stat().st_size / 1e3
    print(f"{len(classes)} accessions, {n_rows} rows -> {args.output} "
          f"({size_kb:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
