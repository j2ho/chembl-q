#!/usr/bin/env python3
"""A ligand is one residue, not every residue sharing its three-letter code.

align_pdb writes every HETATM regardless of chain, so an aligned file holds
each subunit's copy of the ligand. Matching on the code alone pooled them:
Q12791 built a 91-residue "pocket" from 20 copies of a membrane lipid spanning
71.6 A, and 45% of curated_v5 targets were scored on a multi-copy ligand.
"""

import tempfile
from pathlib import Path

import numpy as np

from chembl_curator.protein_filter import ProteinFilter
from chembl_curator.receptor_similarity import chembl_pocket, parse_ligand_atoms

SEPARATION = 60.0     # distance between the two copies' sites


def _atom(serial, name, resname, chain, resseq, xyz, record="ATOM  ", element=None):
    element = element or name[0]
    return (f"{record}{serial:5d} {name:<4s}{resname:>4s}{chain:>2s}{resseq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00"
            f"{element:>12s}\n")


def _site(serial, chain, centre, first_resseq=1):
    """A ligand wrapped in a shell of residues dense enough to count as buried.

    MIN_BURIAL_NEIGHBOURS wants 20 protein heavy atoms within 8 A of the mean
    ligand atom, and pocket_from_structure wants at least 6 residues, so a
    sparse backbone is not enough to exercise either code path.
    """
    centre = np.array(centre, dtype=float)
    lines = []
    resseq = first_resseq
    # 18 residues on two rings around the ligand, 5 A and 7 A out.
    for radius in (5.0, 7.0):
        for k in range(9):
            angle = 2 * np.pi * k / 9
            base = centre + np.array([radius * np.cos(angle),
                                      radius * np.sin(angle),
                                      1.5 * (radius - 6.0)])
            for name, off in (("N", [-1.2, 0, 0]), ("CA", [0, 0, 0]),
                              ("C", [1.2, 0, 0]), ("O", [1.2, 1.0, 0]),
                              ("CB", [0, 0, 1.5])):
                lines.append(_atom(serial, name, "GLY", chain, resseq,
                                   base + np.array(off, dtype=float)))
                serial += 1
            resseq += 1
    return lines, serial, resseq


def _ligand(serial, chain, resseq, centre, n=6):
    centre = np.array(centre, dtype=float)
    lines = []
    for i in range(n):
        xyz = centre + np.array([0.8 * (i - n / 2), 0.0, 0.0])
        lines.append(_atom(serial + i, f"C{i+1}", "LIG", chain, resseq, xyz,
                           record="HETATM", element="C"))
    return lines, serial + n


def _two_copy_structure(path, chain_b="B"):
    """Two sites far apart, each with its own copy of LIG."""
    lines, s, _ = _site(1, "A", [0.0, 0.0, 0.0])
    lig_a, s = _ligand(s, "A", 900, [0.0, 0.0, 0.0])
    lines += lig_a
    more, s, _ = _site(s, chain_b, [0.0, 0.0, SEPARATION])
    lines += more
    lig_b, s = _ligand(s, chain_b, 900, [0.0, 0.0, SEPARATION])
    lines += lig_b
    Path(path).write_text("".join(lines) + "END\n")


def test_named_ligand_alone_pools_every_copy():
    """The old behaviour, kept as the thing the fix has to beat."""
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "s.pdb"
        _two_copy_structure(p)
        pooled = parse_ligand_atoms(p, "LIG")
    assert len(pooled) == 12
    span = np.linalg.norm(pooled[:, None] - pooled[None, :], axis=-1).max()
    assert span > SEPARATION - 1, span


def test_ligand_residue_selects_one_copy():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "s.pdb"
        _two_copy_structure(p)
        one = parse_ligand_atoms(p, "LIG", ligand_residue="A:900")
        other = parse_ligand_atoms(p, "LIG", ligand_residue="B:900")
    assert len(one) == 6 and len(other) == 6
    assert np.linalg.norm(one.mean(0) - other.mean(0)) > SEPARATION - 1


def test_pocket_is_built_from_one_copy_only():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "s.pdb"
        _two_copy_structure(p)
        pooled = chembl_pocket(p, "LIG", pocket_radius=8.0)
        single = chembl_pocket(p, "LIG", pocket_radius=8.0,
                               ligand_residue="A:900")
    assert pooled is not None and single is not None
    # Pooling reaches the far site as well; naming the copy does not.
    assert len(pooled) == 2 * len(single), (len(pooled), len(single))
    assert all(abs(ca[2]) < SEPARATION / 2 for _, ca in single)


def test_ligand_copies_in_one_chain_stay_separate():
    """Two copies in ONE chain must not merge into a molecule between them."""
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "s.pdb"
        lines, s, next_res = _site(1, "A", [0.0, 0.0, 0.0])
        first, s = _ligand(s, "A", 900, [0.0, 0.0, 0.0])
        more, s, _ = _site(s, "A", [0.0, 0.0, SEPARATION], first_resseq=next_res)
        second, s = _ligand(s, "A", 901, [0.0, 0.0, SEPARATION])
        Path(p).write_text("".join(lines + first + more + second) + "END\n")

        pf = ProteinFilter(Path(tmp), log_level="ERROR")
        ligs = pf.get_ligands_from_pdb(p, ["A"])

    keys = sorted(l.residue_key for l in ligs)
    assert keys == ["A:900", "A:901"], keys
    assert all(l.n_heavy == 6 for l in ligs), [l.n_heavy for l in ligs]
    gap = np.linalg.norm(ligs[0].center - ligs[1].center)
    assert gap > SEPARATION - 1, gap


def test_two_sites_are_not_one_pocket():
    """The single-pocket filter must see two sites, not one merged ligand."""
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "s.pdb"
        lines, s, next_res = _site(1, "A", [0.0, 0.0, 0.0])
        first, s = _ligand(s, "A", 900, [0.0, 0.0, 0.0])
        more, s, _ = _site(s, "A", [0.0, 0.0, SEPARATION], first_resseq=next_res)
        second, s = _ligand(s, "A", 901, [0.0, 0.0, SEPARATION])
        Path(p).write_text("".join(lines + first + more + second) + "END\n")

        pf = ProteinFilter(Path(tmp), log_level="ERROR")
        chain_ligands = pf.get_chain_ligand_contacts(p, ["A"])
        assert pf.is_single_ligand_bound(chain_ligands["A"]) is False


def test_first_model_only_in_ligand_parsing():
    """A second model must not double the atoms of a residue."""
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "nmr.pdb"
        m1, s, _ = _site(1, "A", [0.0, 0.0, 0.0])
        lig1, s = _ligand(s, "A", 900, [0.0, 0.0, 0.0])
        # Model 2: same residue numbering, shifted, as an NMR ensemble differs.
        m2, s2, _ = _site(1, "A", [0.0, 0.0, 5.0])
        lig2, s2 = _ligand(s2, "A", 900, [0.0, 0.0, 5.0])
        Path(p).write_text(
            "MODEL        1\n" + "".join(m1 + lig1) + "ENDMDL\n"
            + "MODEL        2\n" + "".join(m2 + lig2) + "ENDMDL\nEND\n")

        pf = ProteinFilter(Path(tmp), log_level="ERROR")
        ligs = pf.get_ligands_from_pdb(p, ["A"])

    assert len(ligs) == 1, [l.residue_key for l in ligs]
    assert ligs[0].n_heavy == 6, ligs[0].n_heavy
    assert abs(ligs[0].center[2]) < 1e-6, ligs[0].center
