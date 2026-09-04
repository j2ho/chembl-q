#!/usr/bin/env python3
"""Chain-specific PDB writing.

P13738 is a homotetramer (A/B/C/D) with LMU sitting at the A/B interface,
so it catches the case where the whole assembly is written instead of the
one chain the pocket was defined on.
"""

import csv
from pathlib import Path
from chembl_curator.protein_filter import ProteinFilter

TARGET = "P13738"


def _chains_in(pdb_file):
    return {line[21:22].strip()
            for line in pdb_file.read_text().splitlines()
            if line.startswith("ATOM") and line[21:22].strip()}


def test_chain_separation():
    curated_dir = Path("test_curated")
    target_dir = curated_dir / TARGET
    target_dir.mkdir(parents=True, exist_ok=True)

    pf = ProteinFilter(curated_dir, log_level="ERROR")
    assert pf.process_target(TARGET), f"{TARGET} was rejected by stage 2"

    aligned_dir = target_dir / "aligned"
    assert aligned_dir.exists(), "no aligned/ directory was produced"

    single = sorted(p for p in aligned_dir.glob("*_?.pdb")
                    if not p.stem.endswith("_aligned"))
    assert single, "no chain-separated files were written"
    for pdb_file in single:
        want = pdb_file.stem.rsplit("_", 1)[1]
        assert _chains_in(pdb_file) == {want}, \
            f"{pdb_file.name} should hold chain {want} alone, holds {_chains_in(pdb_file)}"

    pocket_csv = target_dir / "pocket_info.csv"
    assert pocket_csv.exists(), "pocket_info.csv was not written"
    with pocket_csv.open() as fh:
        rows = list(csv.DictReader(fh))
    assert rows, "pocket_info.csv has no entries"
    assert rows[0]["Ligand_Name"] == "LMU", \
        f"expected LMU as the representative ligand, got {rows[0]['Ligand_Name']}"
    assert rows[0]["Aligned_File"].endswith(f"_{rows[0]['Chain']}.pdb"), \
        "pocket_info points at a file that is not the chain it names"


if __name__ == "__main__":
    test_chain_separation()
    print("ok")
