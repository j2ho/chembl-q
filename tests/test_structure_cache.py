#!/usr/bin/env python3
"""The shared structure cache, and that the local mirror is opt-in."""

import tempfile
from pathlib import Path
from unittest import mock

from chembl_curator.protein_filter import ProteinFilter

PDB_TEXT = (
    "MODRES 1ABC KCX A  156  LYS  carboxylated\n"
    "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 20.00      N\n"
    "END\n"
)


def _pf(tmp, **kw):
    return ProteinFilter(Path(tmp) / "curated", log_level="CRITICAL", **kw)


def test_cache_dir_defaults_beside_the_curated_dir():
    with tempfile.TemporaryDirectory() as tmp:
        pf = _pf(tmp)
        assert pf.cache_dir == Path(tmp) / "pdb_cache"
        assert pf.cache_dir.is_dir()


def test_a_cached_structure_is_copied_instead_of_downloaded():
    with tempfile.TemporaryDirectory() as tmp:
        pf = _pf(tmp)
        (pf.cache_dir / "1abc.pdb").write_text(PDB_TEXT)
        out = Path(tmp) / "out"
        out.mkdir()
        with mock.patch("chembl_curator.protein_filter.requests.get") as get:
            assert pf.download_pdb("1ABC", out)
            get.assert_not_called()
        assert (out / "1abc.pdb").read_text() == PDB_TEXT


def test_the_sidecar_travels_with_the_structure():
    """Without this the modified-residue filter silently stops working."""
    with tempfile.TemporaryDirectory() as tmp:
        pf = _pf(tmp)
        (pf.cache_dir / "9xyz.pdb").write_text(PDB_TEXT)
        (pf.cache_dir / "9xyz.modres").write_text(
            "comp_id\tchain\tseq_id\tparent\nKCX\tA\t156\tLYS\n")
        out = Path(tmp) / "out"
        out.mkdir()
        assert pf.download_pdb("9XYZ", out)
        assert pf.modified_residues(out / "9xyz.pdb") == {"KCX"}


def test_a_download_is_stored_for_the_next_run():
    with tempfile.TemporaryDirectory() as tmp:
        pf = _pf(tmp)
        out = Path(tmp) / "out"
        out.mkdir()
        resp = mock.Mock(text=PDB_TEXT)
        resp.raise_for_status = mock.Mock()
        with mock.patch("chembl_curator.protein_filter.requests.get",
                        return_value=resp):
            assert pf.download_pdb("2DEF", out)
        assert (pf.cache_dir / "2def.pdb").read_text() == PDB_TEXT


def test_an_unavailable_structure_is_not_retried():
    """RCSB holds thousands of these; retrying each run wastes hours."""
    with tempfile.TemporaryDirectory() as tmp:
        pf = _pf(tmp)
        out = Path(tmp) / "out"
        out.mkdir()
        with mock.patch("chembl_curator.protein_filter.requests.get",
                        side_effect=RuntimeError("404")) as get:
            assert not pf.download_pdb("0BAD", out)
            first = get.call_count
        assert (pf.cache_dir / "0bad.miss").exists()
        with mock.patch("chembl_curator.protein_filter.requests.get") as get2:
            assert not pf.download_pdb("0BAD", out)
            get2.assert_not_called()
        assert first > 0


def test_local_mirror_is_off_by_default():
    """pdb_get exists on one cluster and nowhere a released pipeline runs."""
    with tempfile.TemporaryDirectory() as tmp:
        assert _pf(tmp).use_local_mirror is False
        assert _pf(tmp, use_local_mirror=True).use_local_mirror is True


def test_default_path_never_calls_pdb_get():
    with tempfile.TemporaryDirectory() as tmp:
        pf = _pf(tmp)
        out = Path(tmp) / "out"
        out.mkdir()
        resp = mock.Mock(text=PDB_TEXT)
        resp.raise_for_status = mock.Mock()
        with mock.patch("chembl_curator.protein_filter.shutil.which") as which, \
                mock.patch("chembl_curator.protein_filter.requests.get",
                           return_value=resp):
            pf.download_pdb("3GHI", out)
            which.assert_not_called()


if __name__ == "__main__":
    import sys
    import traceback

    failed = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception:
            failed += 1
            print(f"  FAIL  {name}")
            traceback.print_exc()
    print(f"\n{'FAILED' if failed else 'all tests passed'} ({failed} failures)")
    sys.exit(1 if failed else 0)
