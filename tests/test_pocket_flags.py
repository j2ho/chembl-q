"""The cofactor pocket flag, and the distinction the two classes carry.

`unverifiable` means no structure of the target holds a drug-like ligand, so
nothing says where the actives bind and nothing will. `displaced` means one is
present but the cofactor won the per-structure ranking, so the pocket is wrong
now and will move when a release rebuilds it. A reader comparing releases
should expect exactly the `displaced` ones to change.
"""

import pytest

from chembl_curator.pocket_flags import (
    CO_SUBSTRATES,
    MIN_DRUG_LIKE_ATOMS,
    NOT_A_LIGAND,
    pass_pocket_noncofactor,
    pocket_cofactor_status,
    pocket_ligands,
    prosthetic_codes,
)

PROSTHETIC = frozenset({'HEM', 'FAD', 'FMN', 'PLP'})

POCKET_HEADER = ('PDB_ID,Chain,Aligned_File,Ligand_Name,Ligand_Residue,'
                 'Center_X,Center_Y,Center_Z\n')


def write_target(root, uniprot, rows, hetatms=None):
    """rows: (pdb, chain, aligned_file, ligand). hetatms: {file: [(code, n)]}"""
    tdir = root / uniprot
    (tdir / 'aligned').mkdir(parents=True)
    with open(tdir / 'pocket_info.csv', 'w') as handle:
        handle.write(POCKET_HEADER)
        for pdb, chain, aligned, lig in rows:
            handle.write(f'{pdb},{chain},{aligned},{lig},{chain}:1,'
                         f'0.000,0.000,0.000\n')
    for aligned, groups in (hetatms or {}).items():
        lines = []
        for i, (code, n_atoms) in enumerate(groups):
            for j in range(n_atoms):
                lines.append(
                    f'HETATM{j:5d}  C   {code:>3s} A{900 + i:4d}    '
                    f'   0.000   0.000   0.000  1.00  0.00           C\n')
        (tdir / 'aligned' / aligned).write_text(''.join(lines))
    return tdir


def test_ordinary_pocket_is_not_flagged(tmp_path):
    tdir = write_target(tmp_path, 'P1', [('1ABC', 'A', '1abc_A.pdb', 'STI')])
    assert pocket_cofactor_status(tdir, PROSTHETIC) == ''
    assert pass_pocket_noncofactor('') == 1


def test_cofactor_only_with_no_other_ligand_is_unverifiable(tmp_path):
    tdir = write_target(
        tmp_path, 'P2', [('1ABC', 'A', '1abc_A.pdb', 'HEM')],
        {'1abc_A.pdb': [('HEM', 43)]})
    assert pocket_cofactor_status(tdir, PROSTHETIC) == 'unverifiable'
    assert pass_pocket_noncofactor('unverifiable') == 0


def test_cofactor_chosen_over_a_co_bound_drug_is_displaced(tmp_path):
    """The aromatase case: haem won the ranking, exemestane was right there."""
    tdir = write_target(
        tmp_path, 'P3', [('1ABC', 'A', '1abc_A.pdb', 'HEM')],
        {'1abc_A.pdb': [('HEM', 43), ('EXM', 22)]})
    assert pocket_cofactor_status(tdir, PROSTHETIC) == 'displaced'
    assert pass_pocket_noncofactor('displaced') == 0


def test_a_co_substrate_does_not_make_it_displaced(tmp_path):
    """NAD is not a drug site; finding one is no evidence of somewhere better."""
    tdir = write_target(
        tmp_path, 'P4', [('1ABC', 'A', '1abc_A.pdb', 'FAD')],
        {'1abc_A.pdb': [('FAD', 53), ('NAD', 44)]})
    assert pocket_cofactor_status(tdir, PROSTHETIC) == 'unverifiable'


def test_a_detergent_does_not_make_it_displaced(tmp_path):
    """DET split DHODH's chains; a detergent blob is not a druggable site."""
    tdir = write_target(
        tmp_path, 'P5', [('1ABC', 'A', '1abc_A.pdb', 'FMN')],
        {'1abc_A.pdb': [('FMN', 31), ('DET', 20), ('Y01', 30)]})
    assert pocket_cofactor_status(tdir, PROSTHETIC) == 'unverifiable'


def test_a_small_fragment_does_not_make_it_displaced(tmp_path):
    tdir = write_target(
        tmp_path, 'P6', [('1ABC', 'A', '1abc_A.pdb', 'HEM')],
        {'1abc_A.pdb': [('HEM', 43), ('XYZ', MIN_DRUG_LIKE_ATOMS - 1)]})
    assert pocket_cofactor_status(tdir, PROSTHETIC) == 'unverifiable'


def test_one_cofactor_pocket_among_several_does_not_flag(tmp_path):
    """The flag needs EVERY pocket to be a cofactor, not just one of them."""
    tdir = write_target(tmp_path, 'P7', [
        ('1ABC', 'A', '1abc_A.pdb', 'HEM'),
        ('2ABC', 'A', '2abc_A.pdb', 'STI'),
    ])
    assert pocket_cofactor_status(tdir, PROSTHETIC) == ''


def test_drug_in_a_second_structure_still_counts_as_displaced(tmp_path):
    tdir = write_target(
        tmp_path, 'P8',
        [('1ABC', 'A', '1abc_A.pdb', 'HEM'),
         ('2ABC', 'A', '2abc_A.pdb', 'HEM')],
        {'1abc_A.pdb': [('HEM', 43)],
         '2abc_A.pdb': [('HEM', 43), ('EXM', 22)]})
    assert pocket_cofactor_status(tdir, PROSTHETIC) == 'displaced'


def test_pocket_ligands_deduplicates_and_sorts(tmp_path):
    tdir = write_target(tmp_path, 'P9', [
        ('1ABC', 'A', '1abc_A.pdb', 'PMP'),
        ('1ABC', 'B', '1abc_B.pdb', 'PLP'),
        ('2ABC', 'A', '2abc_A.pdb', 'PLP'),
    ])
    assert pocket_ligands(tdir) == ['PLP', 'PMP']


def test_missing_asset_blanks_the_flag_rather_than_passing_it(tmp_path):
    """0 would read as "checked, the pocket is fine"."""
    assert pass_pocket_noncofactor(None) == ''


def test_a_missing_target_directory_is_not_flagged(tmp_path):
    assert pocket_cofactor_status(tmp_path / 'nope', PROSTHETIC) == ''


def test_the_shipped_asset_loads_and_excludes_co_substrates():
    codes = prosthetic_codes()
    assert {'HEM', 'HEC', 'FAD', 'FMN', 'PLP', 'TPP', 'SF4'} <= codes
    assert not (codes & CO_SUBSTRATES)
    assert not (codes & NOT_A_LIGAND)


@pytest.mark.parametrize('status,expected', [
    ('', 1), ('displaced', 0), ('unverifiable', 0)])
def test_polarity_matches_the_other_pass_columns(status, expected):
    """1 is the target clearing the bar, as for pass_reps5."""
    assert pass_pocket_noncofactor(status) == expected
