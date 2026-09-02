#!/usr/bin/env python3
"""Geometry tests for order-free pocket superposition."""

import numpy as np

from chembl_curator.pocket_align import (
    aligned_rmsd,
    fingerprint,
    kabsch,
    pocket_from_structure,
)

RESIDUES = ["ALA", "ARG", "ASP", "SER", "TRP", "LEU", "GLU", "THR"]


def _pocket(seed=0, n=8):
    rng = np.random.default_rng(seed)
    xyz = rng.normal(scale=6.0, size=(n, 3))
    return [(RESIDUES[i % len(RESIDUES)], xyz[i]) for i in range(n)]


def _rotate(pocket, seed=1, shift=(10.0, -4.0, 7.0)):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return [(name, xyz @ q.T + np.array(shift)) for name, xyz in pocket]


def test_rigid_motion_recovers_zero_rmsd():
    p = _pocket()
    n, rmsd = aligned_rmsd(p, _rotate(p))
    assert n >= 6
    assert rmsd < 1e-6, rmsd


def test_shuffled_residue_order_still_matches():
    """The whole point: TM-align cannot do this, Hungarian must."""
    p = _pocket()
    moved = _rotate(p)
    shuffled = [moved[i] for i in [5, 0, 7, 2, 6, 1, 4, 3]]
    n, rmsd = aligned_rmsd(p, shuffled)
    assert rmsd < 1e-6, rmsd


def test_unrelated_pockets_give_large_rmsd():
    n, rmsd = aligned_rmsd(_pocket(seed=0), _pocket(seed=99))
    assert rmsd > 2.0, rmsd


def test_kabsch_rejects_reflections():
    p = _pocket()
    xyz = np.stack([c for _, c in p])
    mirrored = xyz * np.array([1.0, 1.0, -1.0])
    rot, _ = kabsch(mirrored, xyz)
    assert np.linalg.det(rot) > 0


def test_fingerprint_is_rotation_invariant():
    p = _pocket()
    a = np.stack([c for _, c in p])
    b = np.stack([c for _, c in _rotate(p)])
    assert np.allclose(fingerprint(a), fingerprint(b), atol=1e-8)


def test_fingerprint_pads_small_pockets():
    xyz = np.zeros((3, 3))
    xyz[1, 0] = 1.0
    xyz[2, 1] = 2.0
    assert fingerprint(xyz, neighbours=8).shape == (3, 8)


def test_trimming_survives_one_bad_residue():
    """A single displaced residue must not wreck the fit for the rest."""
    p = _pocket(n=12)
    moved = _rotate(p)
    moved[3] = (moved[3][0], moved[3][1] + np.array([40.0, 40.0, 40.0]))
    n, rmsd = aligned_rmsd(p, moved)
    assert n <= 11, "the outlier should have been trimmed"
    assert rmsd < 1e-5, rmsd


def test_pocket_selection_uses_all_heavy_atoms():
    """A residue whose CA is far but whose side chain touches the ligand counts."""
    ligand = np.array([[0.0, 0.0, 0.0]])
    residues = [
        ("ALA", np.array([20.0, 0.0, 0.0]), np.array([[20.0, 0.0, 0.0], [5.0, 0.0, 0.0]])),
        ("SER", np.array([1.0, 0.0, 0.0]), np.array([[1.0, 0.0, 0.0]])),
        ("LEU", np.array([50.0, 0.0, 0.0]), np.array([[50.0, 0.0, 0.0]])),
    ]
    pocket = pocket_from_structure(residues, ligand, radius=8.0, min_residues=1)
    assert [name for name, _ in pocket] == ["ALA", "SER"]


def test_pocket_selection_rejects_tiny_and_empty():
    ligand = np.array([[0.0, 0.0, 0.0]])
    residues = [("SER", np.zeros(3), np.zeros((1, 3)))]
    assert pocket_from_structure(residues, ligand, min_residues=6) is None
    assert pocket_from_structure(residues, np.empty((0, 3))) is None


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
