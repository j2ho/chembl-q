#!/usr/bin/env python3
"""Test PDB download through whichever route is available."""

import tempfile
from pathlib import Path
from chembl_curator import ProteinFilter

RECORD_STARTS = ("HEADER", "ATOM", "HETATM", "REMARK")


def test_pdb_download():
    """1ATP must arrive as a parseable PDB. Requires network."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        pdb_id = "1ATP"

        pf = ProteinFilter(curated_dir=tmpdir, log_level="ERROR")
        assert pf.download_pdb(pdb_id, tmpdir), f"download_pdb failed for {pdb_id}"

        out = tmpdir / f"{pdb_id.lower()}.pdb"
        assert out.exists(), f"download reported success but {out.name} is missing"
        assert out.stat().st_size > 0, f"{out.name} is empty"

        head = out.read_text().splitlines()[:5]
        assert any(line.startswith(RECORD_STARTS) for line in head), \
            f"{out.name} does not begin with any PDB record: {head}"


if __name__ == "__main__":
    test_pdb_download()
    print("ok")
