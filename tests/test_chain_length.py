#!/usr/bin/env python3
"""The receptor chain must be a protein, not a peptide of one.

A UniProt entry often maps to a PDB chain holding only a fragment of it, and
that fragment is usually the ligand of some other protein in the entry. In the
v6 build Q92934 came through as the 6-residue BH3 peptide WAQRGR of a
168-residue protein, and P01106 as 10 residues of 454, so their pockets
described whatever the peptide was bound to. There was an upper bound on chain
size and no lower one, and pocket_from_structure happens to accept 6 residues.

Two bars, because neither alone is enough: the absolute one passes a
23-residue piece of a 1,426-residue protein, and the coverage one alone would
cut a genuine single domain of a large multi-domain protein.
"""

import tempfile
from pathlib import Path

import pytest

from chembl_curator.protein_filter import ProteinFilter


def _chain(chain, n, start=1):
    lines = []
    for i in range(n):
        x = 3.8 * i
        lines.append(
            f"ATOM  {start + i:5d}  CA  GLY{chain:>2s}{i + 1:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n")
    return lines


def _write(path, chains):
    lines, serial = [], 1
    for ch, n in chains:
        lines += _chain(ch, n, serial)
        serial += n
    Path(path).write_text("".join(lines) + "END\n")


@pytest.fixture
def pf(tmp_path):
    return ProteinFilter(tmp_path, log_level="ERROR")


def test_peptide_chain_is_dropped(pf, tmp_path):
    """Q92934's case: 6 residues of a 168-residue protein."""
    p = tmp_path / "s.pdb"
    _write(p, [("M", 6), ("A", 300)])
    assert pf._chains_long_enough(p, ["M"], 168) == []


def test_full_length_chain_is_kept(pf, tmp_path):
    p = tmp_path / "s.pdb"
    _write(p, [("A", 300)])
    assert pf._chains_long_enough(p, ["A"], 330) == ["A"]


def test_fragment_over_the_absolute_bar_still_fails_on_coverage(pf, tmp_path):
    """O00512's case: 23 of 1,426 residues clears nothing, but 60 of 1,426
    clears the absolute bar and is still a fragment."""
    p = tmp_path / "s.pdb"
    _write(p, [("L", 60)])
    assert pf._chains_long_enough(p, ["L"], 1426) == []


def test_single_domain_of_a_large_protein_is_kept(pf, tmp_path):
    """Coverage alone would cut this; the absolute bar is what saves it."""
    p = tmp_path / "s.pdb"
    _write(p, [("A", 196)])
    assert pf._chains_long_enough(p, ["A"], 1202) == ["A"]


def test_only_the_short_chains_are_dropped(pf, tmp_path):
    p = tmp_path / "s.pdb"
    _write(p, [("A", 300), ("P", 8), ("B", 280)])
    assert pf._chains_long_enough(p, ["A", "P", "B"], 320) == ["A", "B"]


def test_counting_is_per_chain_not_summed(pf, tmp_path):
    """Two short chains must not add up to one long enough one."""
    p = tmp_path / "s.pdb"
    _write(p, [("A", 30), ("B", 30)])
    assert pf._chains_long_enough(p, ["A", "B"], 60) == []


def test_later_models_do_not_inflate_the_count(pf, tmp_path):
    """An NMR ensemble would otherwise multiply a peptide past the bar."""
    p = tmp_path / "nmr.pdb"
    body = "".join(_chain("M", 10))
    p.write_text("".join(
        f"MODEL     {i:4d}\n{body}ENDMDL\n" for i in range(1, 9)) + "END\n")
    assert pf._count_residues_per_chain(p) == {"M": 10}
    assert pf._chains_long_enough(p, ["M"], 168) == []


def test_bars_can_be_disabled(tmp_path):
    p = tmp_path / "s.pdb"
    _write(p, [("M", 6)])
    off = ProteinFilter(tmp_path, log_level="ERROR",
                        min_chain_residues=0, min_chain_coverage=0.0)
    assert off._chains_long_enough(p, ["M"], 168) == ["M"]


def test_unknown_uniprot_length_falls_back_to_the_absolute_bar(pf, tmp_path):
    """A missing sequence must not silently disable the coverage check into
    passing everything, nor reject everything."""
    p = tmp_path / "s.pdb"
    _write(p, [("A", 300), ("M", 6)])
    assert pf._chains_long_enough(p, ["A", "M"], 0) == ["A"]
