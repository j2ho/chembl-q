#!/usr/bin/env python3
"""Stage 7: external pocket extraction, caching, and entry resolution."""

import tempfile
from pathlib import Path

import numpy as np

from chembl_curator.external_pockets import (RESIDUE_CODES,
                                             ExternalPocketLeakage,
                                             resolve_external_entry)
from chembl_curator.pocket_align import aligned_rmsd, prepare_pocket
from chembl_curator.receptor_similarity import (parse_ligand_atoms,
                                                parse_protein_residues)


def _write(path, lines):
    path.write_text("".join(lines))
    return path


def _atom(serial, name, resname, chain, resseq, xyz, record="ATOM  ", element=None):
    element = element or name.strip()[0]
    return (f"{record}{serial:5d} {name:<4s} {resname:>3s} {chain}{resseq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          "
            f"{element:>2s}\n")


def _toy_receptor(path, n_residues=8, spacing=2.0):
    """A row of glycines, each with N/CA/C, one per `spacing` angstroms."""
    lines = []
    serial = 1
    for i in range(n_residues):
        base = np.array([i * spacing, 0.0, 0.0])
        for name, offset in (("N", [0, 0.5, 0]), ("CA", [0, 0, 0]), ("C", [0, -0.5, 0])):
            lines.append(_atom(serial, name, "GLY", "A", i + 1, base + offset))
            serial += 1
    return _write(path, lines)


def _toy_ligand(path, xyz=(0.0, 3.0, 0.0)):
    return _write(path, [_atom(1, "C1", "LIG", "A", 1, np.array(xyz),
                               record="HETATM", element="C")])


# ── entry resolution ──────────────────────────────────────────────────────

def test_pdbbind_falls_back_to_the_renamed_ligand():
    """The general set ships only *_ligand_renamed.pdb.

    Requiring the plain name resolved 4,518 of 15,442 PDBbind entries and
    silently left the rest unchecked for leakage.
    """
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "v2020-others"
        entry = root / "13gs"
        entry.mkdir(parents=True)
        (entry / "13gs_protein.pdb").write_text("")
        (entry / "13gs_ligand_renamed.pdb").write_text("")

        got = resolve_external_entry("pdbbind", "13gs", None, [root])
        assert got is not None, "entry with only the renamed ligand did not resolve"
        assert got[1].name == "13gs_ligand_renamed.pdb"


def test_pdbbind_prefers_the_plain_ligand_when_both_exist():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "v2020-refined"
        entry = root / "10gs"
        entry.mkdir(parents=True)
        (entry / "10gs_protein.pdb").write_text("")
        (entry / "10gs_ligand.pdb").write_text("")
        (entry / "10gs_ligand_renamed.pdb").write_text("")

        got = resolve_external_entry("pdbbind", "10gs", None, [root])
        assert got[1].name == "10gs_ligand.pdb"


def test_pdbbind_searches_every_given_set():
    with tempfile.TemporaryDirectory() as tmp:
        refined = Path(tmp) / "refined"
        others = Path(tmp) / "others"
        refined.mkdir()
        entry = others / "1a42"
        entry.mkdir(parents=True)
        (entry / "1a42_protein.pdb").write_text("")
        (entry / "1a42_ligand.pdb").write_text("")

        assert resolve_external_entry("pdbbind", "1a42", None, [refined, others])
        assert resolve_external_entry("pdbbind", "1a42", None, [refined]) is None


def test_biolip_entry_maps_to_its_single_chain_receptor():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "receptor").mkdir()
        (root / "ligand").mkdir()
        (root / "receptor" / "10gsA.pdb").write_text("")
        (root / "ligand" / "10gs_VWW_A_1.pdb").write_text("")

        got = resolve_external_entry("biolip", "10gs_VWW_A_1", root, [])
        assert got is not None
        assert got[0].name == "10gsA.pdb"


def test_missing_files_resolve_to_none_rather_than_a_bad_path():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "receptor").mkdir()
        (root / "ligand").mkdir()
        assert resolve_external_entry("biolip", "9xyz_ABC_A_1", root, []) is None
        assert resolve_external_entry("nonsense", "whatever", root, []) is None


# ── extraction ────────────────────────────────────────────────────────────

def test_pocket_holds_only_residues_near_the_ligand():
    """The ligand sits 3 A off the first residue, so the far end of the chain
    must be excluded at an 8 A radius."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        receptor = _toy_receptor(tmp / "rec.pdb", n_residues=10, spacing=3.0)
        ligand = _toy_ligand(tmp / "lig.pdb", xyz=(0.0, 3.0, 0.0))

        residues = parse_protein_residues(receptor)
        assert len(residues) == 10

        from chembl_curator.pocket_align import pocket_from_structure
        pocket = pocket_from_structure(residues, parse_ligand_atoms(ligand),
                                       radius=8.0, min_residues=1)
        assert pocket is not None
        assert 1 <= len(pocket) < 10, \
            f"expected a local pocket, got {len(pocket)} of 10 residues"


def test_ligand_parsed_without_a_name_takes_every_atom():
    """PDBbind and BioLiP ship the ligand as its own file, and PDBbind writes
    residue name and numbering that no ChEMBL-side name would match."""
    with tempfile.TemporaryDirectory() as tmp:
        ligand = _toy_ligand(Path(tmp) / "lig.pdb")
        assert len(parse_ligand_atoms(ligand)) == 1
        assert len(parse_ligand_atoms(ligand, "LIG")) == 1
        assert len(parse_ligand_atoms(ligand, "XXX")) == 0


# ── cache round trip ──────────────────────────────────────────────────────

def test_cache_round_trip_preserves_the_superposition():
    """A pocket that has been through the npz must give the same answer as one
    that has not, or every stage 7 number is measuring the cache."""
    rng = np.random.default_rng(0)
    names = [RESIDUE_CODES[i % len(RESIDUE_CODES)] for i in range(20)]
    coords = rng.normal(scale=5.0, size=(20, 3))
    pocket = list(zip(names, coords))

    with tempfile.TemporaryDirectory() as tmp:
        cache = Path(tmp) / "ext.npz"
        codes = np.array([RESIDUE_CODES.index(n) for n in names], dtype=np.uint8)
        np.savez_compressed(
            cache,
            keys=np.array(["biolip.test_LIG_A_1"]),
            lengths=np.array([20], dtype=np.int32),
            codes=codes,
            coords=coords.astype(np.float32),
            pocket_radius=np.array(8.0),
        )
        loaded = ExternalPocketLeakage.load(cache)

    assert len(loaded) == 1
    key, prepared = loaded[0]
    assert key == "biolip.test_LIG_A_1"

    # float32 storage is the only difference, so allow its rounding but nothing more.
    n_direct, rmsd_direct = aligned_rmsd(pocket, pocket)
    from chembl_curator.pocket_align import aligned_rmsd_prepared
    n_cached, rmsd_cached = aligned_rmsd_prepared(prepare_pocket(pocket), prepared)
    assert n_direct == n_cached
    assert abs(rmsd_direct - rmsd_cached) < 1e-3, (rmsd_direct, rmsd_cached)


def test_cache_keeps_residue_identity():
    """Residue class drives the mismatch penalty; codes that decode to the
    wrong name would quietly change every RMSD."""
    with tempfile.TemporaryDirectory() as tmp:
        cache = Path(tmp) / "ext.npz"
        names = ["TRP", "GLY", "ASP", "ALA", "MSE", "HIS"]
        np.savez_compressed(
            cache,
            keys=np.array(["pdbbind.1abc"]),
            lengths=np.array([len(names)], dtype=np.int32),
            codes=np.array([RESIDUE_CODES.index(n) for n in names], dtype=np.uint8),
            coords=np.arange(len(names) * 3, dtype=np.float32).reshape(-1, 3),
            pocket_radius=np.array(8.0),
        )
        _, prepared = ExternalPocketLeakage.load(cache)[0]

    from chembl_curator.pocket_align import AA_CLASS
    assert list(prepared[0]) == [AA_CLASS[n] for n in names]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
    print("ok")
