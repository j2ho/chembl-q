"""The cross-active bar must see a target's whole reference set.

Three defects measured on the shipped build, each of which the old test suite
passed over because it only ever exercised cluster representatives:

  * 94 assignments sat above 0.9 similarity to a NON-representative active of
    their own target, one of them an exact duplicate.
  * 14 assignments were the same molecule as a cocrystal ligand of the target
    they were a decoy for. Cocrystal ligands carry no ChEMBL id, so no
    id-based exclusion can see them.
  * weak binders - activity observed above active_max_nm - were in no
    reference set at all, and were not written to a file either.
"""

import gzip
import pickle
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from chembl_curator.compound_pool import CompoundPool
from chembl_curator.decoy_selector import DecoySelector

# Four unrelated scaffolds, pairwise ECFP4 below 0.21 (measured). Each
# reference member can therefore only reject its OWN duplicate, so a passing
# test pins down which reference entry did the work. Picking one scaffold
# family instead made the control vacuous: every candidate was above the bar
# against the representative alone, before any of the new entries mattered.
REPRESENTATIVE = "c1ccc(C(=O)Nc2ccc(Cl)cc2)cc1"
NON_REP_ACTIVE = "Cc1nc(Nc2ccccc2S(N)(=O)=O)nc(N)c1C#N"
WEAK_BINDER = "OC(=O)Cc1c[nH]c2ccccc12"
COCRYSTAL = "CC12CCC3CCCCC3C1CCC2O"
FAR = "CCCCCCCCCCCCCCCC(=O)O"

POOL_SMILES = {
    "REP": REPRESENTATIVE,
    "COPY_OF_NON_REP": NON_REP_ACTIVE,
    "COPY_OF_WEAK": WEAK_BINDER,
    "COPY_OF_XTAL": COCRYSTAL,
    "FARAWAY": FAR,
}

def _fp(smiles):
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(
        Chem.MolFromSmiles(smiles), 2, 2048)


def _entry(smiles, targets):
    mol = Chem.MolFromSmiles(smiles)
    return {
        "smiles": smiles,
        "n_heavy": mol.GetNumHeavyAtoms(),
        "mw": Descriptors.MolWt(mol),
        "clogp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "n_arm_ring": rdMolDescriptors.CalcNumAromaticRings(mol),
        "n_hbd": rdMolDescriptors.CalcNumHBD(mol),
        "n_hba": rdMolDescriptors.CalcNumHBA(mol),
        "fp": _fp(smiles),
        "target_set": targets,
    }


def _build(tmp, reference_smiles):
    """A target whose only representative is REP, with a chosen reference set."""
    data_dir = Path(tmp)
    (data_dir / "P0").mkdir(parents=True)
    (data_dir / "passed_targets.txt").write_text("P0\n")
    pool = {k: _entry(v, {"P0"} if k == "REP" else {"PX"})
            for k, v in POOL_SMILES.items()}
    with open(data_dir / "compound_pool.pkl", "wb") as handle:
        pickle.dump({
            "pool": pool,
            "target_actives": {"P0": {"REP"}},
            "target_all_actives": {"P0": {"REP"}},
            "target_ref_fps": {"P0": [_fp(s) for s in reference_smiles]},
        }, handle)
    return data_dir


def _decoys(data_dir, **kw):
    opts = dict(max_decoys=4, tanimoto_thresh=0.95, cross_active_thresh=0.5,
                mw_window=500.0, clogp_window=10.0, tpsa_window=300.0,
                hbd_window=10, hba_window=10, arm_ring_window=10,
                max_selection_count=99, log_level="ERROR")
    opts.update(kw)
    DecoySelector(**opts).run(data_dir=data_dir)
    rows = (data_dir / "P0" / "decoys.tsv").read_text().splitlines()[1:]
    out = {}
    for row in rows:
        parts = row.split("\t")
        out[parts[0]] = set(parts[1].split(";")) if len(parts) > 1 and parts[1] else set()
    return out


def test_duplicate_of_a_non_representative_active_is_rejected(tmp_path):
    """The 94-assignment defect, in miniature.

    REP is the only representative. NON_REP_ACTIVE is a real active of P0 that
    lost its Butina cluster vote, so the old code never compared against it,
    and COPY_OF_NON_REP is the same molecule under another ChEMBL id.
    """
    data_dir = _build(tmp_path, [REPRESENTATIVE, NON_REP_ACTIVE])
    assert "COPY_OF_NON_REP" not in _decoys(data_dir)["REP"]


def test_representatives_only_reference_lets_it_through(tmp_path):
    """Control: with the old reference set the same candidate is selected.

    Without this the test above could pass for any reason at all.
    """
    data_dir = _build(tmp_path, [REPRESENTATIVE])
    assert "COPY_OF_NON_REP" in _decoys(data_dir)["REP"]


def test_cocrystal_ligand_duplicate_is_rejected(tmp_path):
    """The 14-assignment defect: a decoy identical to a cocrystal ligand."""
    data_dir = _build(tmp_path, [REPRESENTATIVE, COCRYSTAL])
    assert "COPY_OF_XTAL" not in _decoys(data_dir)["REP"]


def test_weak_binder_duplicate_is_rejected(tmp_path):
    """A weak binder is still a binder."""
    data_dir = _build(tmp_path, [REPRESENTATIVE, WEAK_BINDER])
    assert "COPY_OF_WEAK" not in _decoys(data_dir)["REP"]


def test_a_dissimilar_candidate_still_gets_through(tmp_path):
    """The bar must not reject everything; that would pass all of the above."""
    data_dir = _build(tmp_path, [REPRESENTATIVE, NON_REP_ACTIVE, COCRYSTAL,
                                 WEAK_BINDER])
    assert "FARAWAY" in _decoys(data_dir)["REP"]


def test_missing_reference_set_is_refused_not_ignored(tmp_path):
    """An old pool pickle must error rather than silently weaken the bar."""
    data_dir = Path(tmp_path)
    (data_dir / "P0").mkdir(parents=True)
    (data_dir / "passed_targets.txt").write_text("P0\n")
    pool = {k: _entry(v, {"P0"} if k == "REP" else {"PX"})
            for k, v in POOL_SMILES.items()}
    with open(data_dir / "compound_pool.pkl", "wb") as handle:
        pickle.dump({"pool": pool, "target_actives": {"P0": {"REP"}},
                     "target_all_actives": {"P0": {"REP"}}}, handle)
    with pytest.raises(KeyError, match="target_ref_fps"):
        _decoys(data_dir)


# ── The reference set builder ────────────────────────────────────────────────

def _pocket_info(rows):
    head = ("PDB_ID,Chain,Aligned_File,Ligand_Name,Ligand_Residue,"
            "Center_X,Center_Y,Center_Z\n")
    return head + "".join(
        f"{pdb.upper()},A,{pdb}_A.pdb,{lig},A:{n},{x},{y},{z}\n"
        for n, (pdb, lig, x, y, z) in enumerate(rows, 1)
    )


def _target(tmp, rows, best="1abc_A"):
    data_dir = Path(tmp)
    tdir = data_dir / "P0"
    tdir.mkdir(parents=True)
    (data_dir / "best_structure.tsv").write_text(
        f"uniprot\tpdbid_chain\tresolution\nP0\t{best}\t2.00\n")
    (tdir / "pocket_info.csv").write_text(_pocket_info(rows))
    return tdir


def _real_comp_ids(n):
    """Component ids that actually resolve, taken from the bundled asset."""
    asset = (Path(__file__).parent.parent / "chembl_curator" / "assets" /
             "ccd_smiles.tsv.gz")
    if not asset.exists():
        pytest.skip("ccd_smiles.tsv.gz not built")
    out = []
    with gzip.open(asset, "rt") as handle:
        next(handle)
        for line in handle:
            comp_id, smiles = line.rstrip("\n").split("\t")[:2]
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None and mol.GetNumHeavyAtoms() > 12:
                out.append(comp_id)
            if len(out) == n:
                break
    return out


def test_every_pdb_contributes_its_pocket_ligand(tmp_path):
    """Several PDB entries of one target, each with a ligand in the pocket.

    All of them are solved structures of something bound here, so all of them
    belong in the reference set, not just the representative's.
    """
    a, b, c = _real_comp_ids(3)
    tdir = _target(tmp_path, [
        ("1abc", a, 0.0, 0.0, 0.0),
        ("2abc", b, 1.0, 0.0, 0.0),
        ("3abc", c, 0.0, 1.5, 0.0),
    ])
    assert len(CompoundPool(log_level="ERROR")._cocrystal_fps(tdir)) == 3


def test_oligomer_copies_collapse_to_one(tmp_path):
    """One ligand repeated per chain must not be counted several times."""
    a, = _real_comp_ids(1)
    tdir = _target(tmp_path, [
        ("1abc", a, 0.0, 0.0, 0.0),
        ("1abc", a, 1.0, 1.0, 1.0),
        ("1abc", a, 2.0, 0.0, 1.0),
    ])
    assert len(CompoundPool(log_level="ERROR")._cocrystal_fps(tdir)) == 1


def test_a_second_ligand_elsewhere_in_one_pdb_is_dropped(tmp_path):
    """The case to be careful about: one PDB, several ligands, one off-site.

    The distance guard is measured from the representative structure's ligand,
    so an off-site ligand in that same entry cannot drag the reference point
    onto itself.
    """
    pocket, far = _real_comp_ids(2)
    tdir = _target(tmp_path, [
        ("1abc", pocket, 0.0, 0.0, 0.0),
        ("1abc", far, 40.0, 0.0, 0.0),
    ])
    assert len(CompoundPool(log_level="ERROR")._cocrystal_fps(tdir)) == 1


def test_a_ligand_in_the_pocket_in_any_copy_is_kept(tmp_path):
    """Off-site in one entry, in the pocket in another: keep it."""
    pocket, other = _real_comp_ids(2)
    tdir = _target(tmp_path, [
        ("1abc", pocket, 0.0, 0.0, 0.0),
        ("2abc", other, 40.0, 0.0, 0.0),
        ("3abc", other, 1.0, 0.0, 0.0),
    ])
    assert len(CompoundPool(log_level="ERROR")._cocrystal_fps(tdir)) == 2


# ── The receptor-similarity exclusion ────────────────────────────────────────

def test_non_representative_active_of_a_similar_receptor_is_excluded(tmp_path):
    """The 13,030-assignment leak.

    SIMPAL is a documented active of PX, a receptor similar enough to P0 to be
    excluded, but it is not one of PX's cluster representatives. It sits in the
    pool because it represents a cluster of some third target, exactly the
    9,479-compound overlap the measurement found. Nothing else catches it: the
    cross-active bar compares against P0's reference set, not PX's.
    """
    data_dir = Path(tmp_path)
    (data_dir / "P0").mkdir(parents=True)
    (data_dir / "passed_targets.txt").write_text("P0\nPX\n")
    (data_dir / "PX").mkdir()
    # P0 and PX are similar receptors.
    (data_dir / "pairwise_seqid.tsv").write_text(
        "query\ttarget\tseqid\nP0\tPX\t0.90\n")

    pool = {k: _entry(v, {"P0"} if k == "REP" else {"PT"})
            for k, v in POOL_SMILES.items()}
    with open(data_dir / "compound_pool.pkl", "wb") as handle:
        pickle.dump({
            "pool": pool,
            # PX's representative is something else entirely; FARAWAY is one of
            # its actives without being its representative.
            "target_actives": {"P0": {"REP"}, "PX": {"COPY_OF_XTAL"}},
            "target_all_actives": {"P0": {"REP"},
                                   "PX": {"COPY_OF_XTAL", "FARAWAY"}},
            "target_ref_fps": {"P0": [_fp(REPRESENTATIVE)], "PX": []},
        }, handle)

    selected = _decoys(data_dir, seqid_thresh=0.6)["REP"]
    assert "FARAWAY" not in selected


def test_representatives_only_exclusion_lets_it_through(tmp_path):
    """Control: keyed on representatives, the same candidate is selected."""
    data_dir = Path(tmp_path)
    (data_dir / "P0").mkdir(parents=True)
    (data_dir / "passed_targets.txt").write_text("P0\nPX\n")
    (data_dir / "PX").mkdir()
    (data_dir / "pairwise_seqid.tsv").write_text(
        "query\ttarget\tseqid\nP0\tPX\t0.90\n")

    pool = {k: _entry(v, {"P0"} if k == "REP" else {"PT"})
            for k, v in POOL_SMILES.items()}
    with open(data_dir / "compound_pool.pkl", "wb") as handle:
        pickle.dump({
            "pool": pool,
            "target_actives": {"P0": {"REP"}, "PX": {"COPY_OF_XTAL"}},
            # The old behaviour: only PX's representatives are known.
            "target_all_actives": {"P0": {"REP"}, "PX": {"COPY_OF_XTAL"}},
            "target_ref_fps": {"P0": [_fp(REPRESENTATIVE)], "PX": []},
        }, handle)

    assert "FARAWAY" in _decoys(data_dir, seqid_thresh=0.6)["REP"]


def test_missing_all_actives_is_refused_not_ignored(tmp_path):
    data_dir = _build(tmp_path, [REPRESENTATIVE])
    with open(data_dir / "compound_pool.pkl", "rb") as handle:
        blob = pickle.load(handle)
    del blob["target_all_actives"]
    with open(data_dir / "compound_pool.pkl", "wb") as handle:
        pickle.dump(blob, handle)
    with pytest.raises(KeyError, match="target_all_actives"):
        _decoys(data_dir)
