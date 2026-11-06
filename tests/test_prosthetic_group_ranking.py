"""A prosthetic group must not define a pocket when a real ligand is there.

Contact count and atom count both favour a buried cofactor, so the previous
ranking picked haem over exemestane on aromatase and FAD over imatinib on
NQO2. The pocket then described the cofactor site, a median 8.5 A from the
ligand it beat, against an 8 A pocket radius.
"""

import pytest

from chembl_curator.protein_filter import (
    PROSTHETIC_GROUPS,
    select_representative_ligand,
)


class FakeLigand:
    def __init__(self, name, n_contacts, n_heavy):
        self.ligand_name = name
        self.n_contacts = n_contacts
        self.n_heavy = n_heavy

    def __repr__(self):
        return f"{self.ligand_name}(c={self.n_contacts},h={self.n_heavy})"


# The real numbers' shape: the cofactor wins on both axes.
HAEM = FakeLigand("HEM", n_contacts=120, n_heavy=43)
FAD = FakeLigand("FAD", n_contacts=150, n_heavy=53)
INHIBITOR = FakeLigand("EXM", n_contacts=40, n_heavy=22)
SECOND_INHIBITOR = FakeLigand("STI", n_contacts=55, n_heavy=37)


def test_inhibitor_beats_haem_despite_fewer_contacts_and_atoms():
    assert select_representative_ligand([HAEM, INHIBITOR]) is INHIBITOR


def test_inhibitor_beats_fad():
    assert select_representative_ligand([FAD, INHIBITOR]) is INHIBITOR


def test_order_does_not_decide_it():
    assert select_representative_ligand([INHIBITOR, HAEM]) is INHIBITOR
    assert select_representative_ligand([HAEM, INHIBITOR]) is INHIBITOR


def test_among_two_real_ligands_contacts_still_decide():
    """Demotion must not disturb the ranking it was not meant to change."""
    assert (select_representative_ligand([INHIBITOR, SECOND_INHIBITOR])
            is SECOND_INHIBITOR)


def test_cofactor_alone_still_defines_the_pocket():
    """Something has to define it; 27 targets have nothing else."""
    assert select_representative_ligand([HAEM]) is HAEM
    assert select_representative_ligand([HAEM, FAD]) is FAD


def test_co_substrates_are_not_demoted():
    """NAD, SAM and ATP turn over, and 1T2F wants NAD to win.

    Whether an inhibitor competes with a co-substrate or binds beside it is
    target-specific, so demoting them would move those pockets for no
    measured reason.
    """
    for code in ("NAD", "NAP", "NDP", "SAM", "SAH", "ATP", "ADP", "COA"):
        assert code not in PROSTHETIC_GROUPS
    nad = FakeLigand("NAD", n_contacts=90, n_heavy=44)
    oxamate = FakeLigand("OXM", n_contacts=30, n_heavy=6)
    assert select_representative_ligand([nad, oxamate]) is nad


@pytest.mark.parametrize("code", ["HEM", "HEC", "FAD", "FMN", "PLP", "TPP",
                                  "SF4", "FES"])
def test_the_classes_we_measured_are_all_listed(code):
    assert code in PROSTHETIC_GROUPS
