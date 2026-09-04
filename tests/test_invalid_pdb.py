#!/usr/bin/env python3
"""Test error handling for invalid PDB ID."""

import tempfile
from pathlib import Path
from chembl_curator import ProteinFilter


def test_invalid_pdb():
    """An unknown PDB ID must be reported as a failure, not half-written."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        invalid_pdb_id = "XXXX"

        pf = ProteinFilter(curated_dir=tmpdir, log_level="ERROR")
        success = pf.download_pdb(invalid_pdb_id, tmpdir)

        assert not success, "download_pdb reported success for a nonexistent entry"
        assert not (tmpdir / f"{invalid_pdb_id.lower()}.pdb").exists(), \
            "a file was left behind for a nonexistent entry"


def test_invalid_pdb_leaves_a_miss_marker():
    """The cache must remember the miss, or every run re-asks RCSB for it."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        cache = tmpdir / "cache"
        cache.mkdir()

        pf = ProteinFilter(curated_dir=tmpdir, log_level="ERROR", cache_dir=cache)
        assert not pf.download_pdb("XXXX", tmpdir)
        assert (cache / "xxxx.miss").exists(), \
            "no .miss marker written; the failed lookup will be repeated forever"


if __name__ == "__main__":
    test_invalid_pdb()
    test_invalid_pdb_leaves_a_miss_marker()
    print("ok")
