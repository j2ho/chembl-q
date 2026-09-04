#!/usr/bin/env python3
"""Test the RCSB web route, which is the only one a released install has."""

import tempfile
from pathlib import Path
from unittest.mock import patch
from chembl_curator import ProteinFilter

RECORD_STARTS = ("HEADER", "ATOM", "HETATM", "REMARK")


def test_rcsb_download():
    """With the site-local mirror hidden, the network path must still work.

    This is the route every user outside this machine takes, so it is worth
    exercising even though pdb_get happens to be installed here.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        pdb_id = "1ATP"

        pf = ProteinFilter(curated_dir=tmpdir, log_level="ERROR")
        with patch("shutil.which", return_value=None):
            success = pf.download_pdb(pdb_id, tmpdir)

        assert success, "RCSB fallback failed with pdb_get unavailable"

        out = tmpdir / f"{pdb_id.lower()}.pdb"
        assert out.exists(), f"download reported success but {out.name} is missing"

        head = out.read_text().splitlines()[:5]
        assert any(line.startswith(RECORD_STARTS) for line in head), \
            f"{out.name} does not begin with any PDB record: {head}"


if __name__ == "__main__":
    test_rcsb_download()
    print("ok")
