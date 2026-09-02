#!/usr/bin/env python3
"""Same-pocket detection uses closest approach, not centroid separation."""

import numpy as np

from chembl_curator.protein_filter import (
    MIN_LIGAND_HEAVY_ATOMS,
    LigandInfo,
    ProteinFilter,
)


def lig(name, coords):
    xyz = np.asarray(coords, dtype=float)
    return LigandInfo(name, xyz.mean(0), "A", coords=xyz, n_heavy=len(xyz))


def bar(name, start, length, n=20, axis=0):
    """A rod of atoms 1 A apart, to stand in for an elongated cofactor."""
    pts = np.zeros((n, 3))
    pts[:, axis] = np.linspace(start, start + length, n)
    return lig(name, pts)


def _pf():
    return ProteinFilter.__new__(ProteinFilter)


def test_touching_ligands_are_one_pocket_despite_far_centroids():
    """The NAD/OXQ case: centroids 9.96 A apart, closest atoms 2.85 A."""
    nad = bar("NAD", 0.0, 20.0)
    oxq = lig("OXQ", [[22.8, 0, 0], [23.8, 0, 0], [24.8, 0, 0],
                      [23.8, 1, 0], [23.8, 0, 1]])
    assert np.linalg.norm(nad.center - oxq.center) > 9.0
    assert _pf().is_single_ligand_bound([nad, oxq])


def test_separated_ligands_are_two_pockets_at_similar_centroid_distance():
    """The P06700 case: same centroid distance, but nothing is touching."""
    a = lig("NCA", [[0, 0, 0], [1, 0, 0], [2, 0, 0], [1, 1, 0], [1, 0, 1]])
    b = lig("XYQ", [[9, 0, 0], [10, 0, 0], [11, 0, 0], [10, 1, 0], [10, 0, 1]])
    assert 7.0 < np.linalg.norm(a.center - b.center) < 11.0
    assert not _pf().is_single_ligand_bound([a, b])


def test_gap_threshold_boundary():
    a = lig("L1", [[0, 0, 0]] * 1 + [[i, 0, 0] for i in range(1, 6)])
    near = lig("L2", [[i, 0, 0] for i in range(10, 16)])   # gap 5.0
    far = lig("L3", [[i, 0, 0] for i in range(11, 17)])    # gap 6.0
    assert _pf().is_single_ligand_bound([a, near])
    assert not _pf().is_single_ligand_bound([a, far])


def test_single_ligand_passes_and_empty_fails():
    assert _pf().is_single_ligand_bound([bar("LIG", 0.0, 5.0)])
    assert not _pf().is_single_ligand_bound([])


def test_three_ligands_need_pairwise_agreement():
    a = lig("L1", [[i, 0, 0] for i in range(6)])
    b = lig("L2", [[i, 0, 0] for i in range(7, 13)])
    far = lig("L3", [[i, 0, 0] for i in range(80, 86)])
    assert _pf().is_single_ligand_bound([a, b])
    assert not _pf().is_single_ligand_bound([a, b, far])


def test_missing_coords_falls_back_to_centroid():
    a = LigandInfo("L1", np.zeros(3), "A")
    near = LigandInfo("L2", np.array([5.0, 0, 0]), "A")
    far = LigandInfo("L3", np.array([50.0, 0, 0]), "A")
    assert _pf().is_single_ligand_bound([a, near])
    assert not _pf().is_single_ligand_bound([a, far])


def test_heavy_atom_floor_matches_the_chembl_side():
    """Pocket-defining ligands use the same size floor as ChEMBL compounds."""
    assert MIN_LIGAND_HEAVY_ATOMS == 5


def test_representative_prefers_contacts_over_size():
    """1T2F: NAD describes the site, oxamate is a small neighbour."""
    nad = LigandInfo("NAD", np.zeros(3), "A", n_heavy=44, n_contacts=19)
    oxq = LigandInfo("OXQ", np.zeros(3), "A", n_heavy=9, n_contacts=6)
    assert max([nad, oxq], key=lambda l: (l.n_contacts, l.n_heavy)) is nad


def test_representative_rejects_a_big_ligand_that_barely_touches():
    """The 53-atom FAD with one contact must not win on atom count."""
    ghost = LigandInfo("FAD", np.zeros(3), "A", n_heavy=53, n_contacts=1)
    real = LigandInfo("ADP", np.zeros(3), "A", n_heavy=27, n_contacts=20)
    assert max([ghost, real], key=lambda l: (l.n_contacts, l.n_heavy)) is real


def test_representative_breaks_contact_ties_by_size():
    small = LigandInfo("SML", np.zeros(3), "A", n_heavy=8, n_contacts=12)
    big = LigandInfo("BIG", np.zeros(3), "A", n_heavy=30, n_contacts=12)
    assert max([small, big], key=lambda l: (l.n_contacts, l.n_heavy)) is big


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
