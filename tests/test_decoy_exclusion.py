#!/usr/bin/env python3
"""Decoy selection must never pick a compound already measured against the target."""

import pickle
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from chembl_curator.decoy_selector import DecoySelector

# Chemically unrelated but similar in size, so property matching cannot be the
# thing that rejects them. Tanimoto between them stays well under 0.3.
SMILES = {
    "ACT1": "CCOc1ccc(cc1)C(=O)Nc1ccccc1",
    "MEAS1": "CN1CCN(CC1)c1ncccn1",
    "MEAS2": "OCC1OC(O)C(O)C(O)C1O",
    "FREE1": "CCCCCCCCCC(=O)NCC",
    "FREE2": "CC(C)CCCNC(=O)CCCC",
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


def _build(tmp, write_measured):
    data_dir = Path(tmp)
    (data_dir / "P0").mkdir(parents=True)
    (data_dir / "passed_targets.txt").write_text("P0\n")
    if write_measured:
        (data_dir / "P0" / "measured.tsv").write_text("chembl_id\nACT1\nMEAS1\nMEAS2\n")
    pool = {k: _entry(v, {"P0"} if k == "ACT1" else {"PX"}) for k, v in SMILES.items()}
    with open(data_dir / "compound_pool.pkl", "wb") as f:
        pickle.dump({"pool": pool, "target_actives": {"P0": {"ACT1"}}}, f)
    return data_dir


def _decoys_for(data_dir):
    DecoySelector(max_decoys=4, tanimoto_thresh=0.5, mw_window=500.0,
                  clogp_window=10.0, tpsa_window=300.0, hbd_window=10,
                  hba_window=10, arm_ring_window=10, max_selection_count=99,
                  log_level="ERROR").run(data_dir)
    rows = (data_dir / "P0" / "decoys.tsv").read_text().splitlines()[1:]
    parts = rows[0].split("\t")
    return set(parts[1].split(";")) if len(parts) > 1 and parts[1] else set()


def test_measured_compounds_are_never_decoys():
    with tempfile.TemporaryDirectory() as tmp:
        picked = _decoys_for(_build(tmp, write_measured=True))
    assert "MEAS1" not in picked and "MEAS2" not in picked, picked
    assert "ACT1" not in picked
    assert picked, "property windows were too tight; nothing was selectable"


def test_without_measured_file_those_compounds_are_selectable():
    """Guards the test above: it must be measured.tsv doing the work."""
    with tempfile.TemporaryDirectory() as tmp:
        picked = _decoys_for(_build(tmp, write_measured=False))
    assert {"MEAS1", "MEAS2"} & picked, picked


def test_load_measured_reports_missing_files():
    with tempfile.TemporaryDirectory() as tmp:
        data_dir = _build(tmp, write_measured=True)
        sel = DecoySelector(log_level="ERROR")
        assert sel._load_measured(data_dir, ["P0"]) == {"P0": {"ACT1", "MEAS1", "MEAS2"}}
        assert sel._load_measured(data_dir, ["NOPE"]) == {}


POCKET_TSV = (
    "target_a\ttarget_b\ttm_score\tpocket_rmsd\tn_matched\tn_pocket_a\tn_pocket_b\n"
    "P0\tP1\t-1.0\t1.20\t40\t45\t44\n"   # close fit, many residues -> similar
    "P0\tP2\t-1.0\t1.20\t8\t9\t10\n"     # same RMSD, too few residues -> not
    "P0\tP3\t-1.0\t3.50\t40\t45\t44\n"   # enough residues, RMSD too high -> not
    "P0\tP4\t-1.0\t-1.0\t0\t45\t44\n"    # alignment failed -> not
)


def _pocket_file(tmp, text=POCKET_TSV):
    p = Path(tmp) / "pairwise_pocket.tsv"
    p.write_text(text)
    return p


def test_pocket_similarity_needs_both_rmsd_and_residue_count():
    with tempfile.TemporaryDirectory() as tmp:
        sel = DecoySelector(pocket_rmsd_thresh=2.0, min_matched_residues=15,
                            log_level="ERROR")
        sim = sel._load_pocket_similar(_pocket_file(tmp))
    assert sim.get("P0") == {"P1"}, sim
    assert "P2" not in sim and "P3" not in sim and "P4" not in sim


def test_residue_count_filter_can_be_disabled():
    """Guards the test above: with the floor at 0 the small pocket comes back."""
    with tempfile.TemporaryDirectory() as tmp:
        sel = DecoySelector(pocket_rmsd_thresh=2.0, min_matched_residues=0,
                            log_level="ERROR")
        sim = sel._load_pocket_similar(_pocket_file(tmp))
    assert sim.get("P0") == {"P1", "P2"}, sim


def test_legacy_tsv_without_n_matched_still_loads():
    """Old stage-5 output has 5 columns; it must not crash, only warn."""
    legacy = (
        "target_a\ttarget_b\ttm_score\tpocket_rmsd\n"
        "P0\tP1\t0.9\t1.20\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        sel = DecoySelector(pocket_rmsd_thresh=2.0, min_matched_residues=15,
                            log_level="CRITICAL")
        sim = sel._load_pocket_similar(_pocket_file(tmp, legacy))
    assert sim.get("P0") == {"P1"}


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
