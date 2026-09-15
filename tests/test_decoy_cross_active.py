#!/usr/bin/env python3
"""A decoy must not be a known binder of the target it is a decoy for.

The similarity filter only ever compared a candidate against the one active it
was being paired with, so a decoy for active A was free to be the same
molecule as active B of the same target. measured.tsv catches that only when
the compound happens to have been tested against the target.

The two bars are deliberately different. Against the paired active the bar is
tanimoto_thresh (0.3). Against every other active it is cross_active_thresh
(0.9), because at 0.3 the cut removes 13.2% of what the pool can supply and
refilling needs twice the spare capacity the reuse cap leaves.
"""

import pickle
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from chembl_curator.decoy_selector import DecoySelector

# ACT_S and TRAP are enantiomers: ECFP4 is built without useChirality, so they
# are indistinguishable to the filter even though their InChIKeys differ. This
# is the class that survived every other check in curated_v5.
SMILES = {
    "ACT_A": "CCOc1ccc(cc1)C(=O)Nc1ccccc1",
    "ACT_S": "C[C@H](N)C(=O)Nc1ccc(Cl)cc1",
    "TRAP":  "C[C@@H](N)C(=O)Nc1ccc(Cl)cc1",
    "FREE1": "CCCCCCCCCC(=O)NCC",
    "FREE2": "CC(C)CCCNC(=O)CCCC",
    "FREE3": "CCCCOC(=O)CCCCC",
}


def _entry(smiles, targets):
    mol = Chem.MolFromSmiles(smiles)
    return {
        "smiles": smiles,
        "n_heavy": mol.GetNumHeavyAtoms(),
        "mw": Descriptors.MolWt(mol),
        "clogp": Crippen.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "n_arm_ring": rdMolDescriptors.CalcNumAromaticRings(mol),
        "n_hbd": rdMolDescriptors.CalcNumHBD(mol),
        "n_hba": rdMolDescriptors.CalcNumHBA(mol),
        "fp": rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
        "target_set": targets,
    }


def _build(tmp):
    """P0 has two actives; TRAP is the enantiomer of the second one.

    TRAP carries no ChEMBL record against P0, so neither the ID exclusion nor
    measured.tsv can stop it being drawn as a decoy.
    """
    data_dir = Path(tmp)
    (data_dir / "P0").mkdir(parents=True)
    (data_dir / "passed_targets.txt").write_text("P0\n")
    actives = {"ACT_A", "ACT_S"}
    pool = {k: _entry(v, {"P0"} if k in actives else {"PX"})
            for k, v in SMILES.items()}
    with open(data_dir / "compound_pool.pkl", "wb") as f:
        pickle.dump({"pool": pool, "target_actives": {"P0": actives}}, f)
    return data_dir


def _decoys(data_dir, **kw):
    opts = dict(max_decoys=4, tanimoto_thresh=0.3, mw_window=500.0,
                clogp_window=10.0, tpsa_window=300.0, hbd_window=10,
                hba_window=10, arm_ring_window=10, max_selection_count=99,
                log_level="ERROR")
    opts.update(kw)
    DecoySelector(**opts).run(data_dir)
    out = {}
    for row in (data_dir / "P0" / "decoys.tsv").read_text().splitlines()[1:]:
        parts = row.split("\t")
        out[parts[0]] = (parts[1].split(";")
                         if len(parts) > 1 and parts[1] else [])
    return out


def test_enantiomer_of_another_active_is_not_a_decoy():
    with tempfile.TemporaryDirectory() as tmp:
        picked = _decoys(_build(tmp))
    chosen = {c for ids in picked.values() for c in ids}
    assert "TRAP" not in chosen, picked
    assert chosen, "nothing was selected, so the test proves nothing"


def test_a_loose_cross_active_bar_lets_it_back_in():
    """Confirms the guard is the cross-active bar, not some other filter."""
    with tempfile.TemporaryDirectory() as tmp:
        picked = _decoys(_build(tmp), cross_active_thresh=1.01)
    chosen = {c for ids in picked.values() for c in ids}
    assert "TRAP" in chosen, picked


def test_paired_active_bar_stays_strict():
    """The 0.3 bar against the paired active must still apply."""
    with tempfile.TemporaryDirectory() as tmp:
        data_dir = _build(tmp)
        picked = _decoys(data_dir, cross_active_thresh=1.01,
                         tanimoto_thresh=0.0)
    # Nothing can clear a 0.0 bar against its own active.
    assert all(not ids for ids in picked.values()), picked
