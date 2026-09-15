#!/usr/bin/env python3
"""Stage 6 must not depend on set iteration order, and must honour an explicit cap.

Seeding random covers the pool shuffle but not the order actives are processed
in: that came from list()-ing a set of ID strings, whose order is salted by
PYTHONHASHSEED. Order matters because sel_count accumulates, so whoever draws
first draws from a fuller pool.

The fixture is built so the cap binds hard (32 slots wanted, 16 available).
Under a sorted order the same four actives always win; under hash order the
winners change with the seed, which is what these tests catch.
"""

import pickle
import subprocess
import sys
import tempfile
import textwrap
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from chembl_curator.decoy_selector import DecoySelector

ACTIVES = {
    "CHEMBL_A3": "CCOc1ccc(cc1)C(=O)Nc1ccccc1",
    "CHEMBL_A1": "CN1CCN(CC1)c1ncccn1",
    "CHEMBL_A2": "OCC1OC(O)C(O)C(O)C1O",
    "CHEMBL_A4": "Clc1ccc(cc1)S(=O)(=O)N1CCOCC1",
    "CHEMBL_A5": "COc1cccc(c1)C(=O)NCC1CC1",
    "CHEMBL_A6": "CCN(CC)C(=O)c1ccc(F)cc1",
    "CHEMBL_A7": "OC(=O)Cc1cccc2ccccc12",
    "CHEMBL_A8": "CC(C)Nc1nc(N)nc(n1)N",
}
POOL_EXTRA = {
    "CHEMBL_D1": "CCCCCCCCCC(=O)NCC",
    "CHEMBL_D2": "CC(C)CCCNC(=O)CCCC",
    "CHEMBL_D3": "CCCCOC(=O)CCCCC",
    "CHEMBL_D4": "CCCCCCNC(=O)OCC",
    "CHEMBL_D5": "CCCCCCCCOC(=O)CC",
    "CHEMBL_D6": "CCCCN(CC)C(=O)CCC",
    "CHEMBL_D7": "CCOCCOCCOCCO",
    "CHEMBL_D8": "CCCCCCCCCCO",
}

# 8 actives x 4 decoys = 32 slots wanted; 8 candidates x cap 2 = 16 available.
CAP = 2
MAX_DECOYS = 4


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
    data_dir = Path(tmp)
    (data_dir / "P0").mkdir(parents=True)
    (data_dir / "passed_targets.txt").write_text("P0\n")
    pool = {k: _entry(v, {"P0"}) for k, v in ACTIVES.items()}
    pool.update({k: _entry(v, {"PX"}) for k, v in POOL_EXTRA.items()})
    with open(data_dir / "compound_pool.pkl", "wb") as f:
        pickle.dump({"pool": pool, "target_actives": {"P0": set(ACTIVES)}}, f)
    return data_dir


# Property and similarity filters are opened up so the cap and the exclusion
# set are the only things that can reject a candidate.
_OPTS = dict(max_decoys=MAX_DECOYS, tanimoto_thresh=0.95, mw_window=1000.0,
             clogp_window=100.0, tpsa_window=1000.0, hbd_window=100,
             hba_window=100, arm_ring_window=100, max_selection_count=CAP,
             log_level="ERROR")


def _parse(data_dir):
    rows = (data_dir / "P0" / "decoys.tsv").read_text().splitlines()[1:]
    out = []
    for row in rows:
        parts = row.split("\t")
        out.append((parts[0],
                    parts[1].split(";") if len(parts) > 1 and parts[1] else []))
    return out


def _run(data_dir, **kw):
    opts = dict(_OPTS)
    opts.update(kw)
    DecoySelector(**opts).run(data_dir)
    return _parse(data_dir)


def test_actives_are_processed_in_sorted_order():
    with tempfile.TemporaryDirectory() as tmp:
        result = _run(_build(tmp))
    assert [a for a, _ in result] == sorted(ACTIVES)


def test_output_is_identical_across_hash_seeds():
    """The real regression: same inputs, different PYTHONHASHSEED, same file.

    Run in subprocesses because PYTHONHASHSEED is read at interpreter start
    and cannot be changed from inside a running one.
    """
    driver = textwrap.dedent(f"""
        import sys
        from pathlib import Path
        sys.path.insert(0, {str(Path(__file__).resolve().parent.parent)!r})
        from chembl_curator.decoy_selector import DecoySelector
        DecoySelector(**{_OPTS!r}).run(Path(sys.argv[1]))
        sys.stdout.write((Path(sys.argv[1]) / "P0" / "decoys.tsv").read_text())
    """)
    outputs = []
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "driver.py"
        script.write_text(driver)
        for seed in ("0", "1", "2", "3"):
            run_dir = Path(tmp) / f"run{seed}"
            _build(run_dir)
            proc = subprocess.run(
                [sys.executable, str(script), str(run_dir)],
                capture_output=True, text=True,
                env={"PYTHONHASHSEED": seed, "PATH": "/usr/bin:/bin"},
            )
            assert proc.returncode == 0, proc.stderr
            outputs.append(proc.stdout)

    assert len(set(outputs)) == 1, (
        "decoy assignment changed with PYTHONHASHSEED:\n"
        + "\n---\n".join(outputs))
    # The fixture only exercises ordering if the cap actually bit.
    assert outputs[0].count(";") > 0


def test_repeated_runs_agree():
    with tempfile.TemporaryDirectory() as tmp:
        data_dir = _build(tmp)
        first = _run(data_dir)
        second = _run(data_dir)
    assert first == second


def test_explicit_cap_is_honoured():
    """An explicit max_selection_count must bind, not be recomputed away."""
    with tempfile.TemporaryDirectory() as tmp:
        result = _run(_build(tmp))
    used = Counter(d for _, decoys in result for d in decoys)
    assert used, "no decoys were selected, so the cap was not exercised"
    assert max(used.values()) <= CAP, dict(used)


def test_cap_below_demand_underfills_rather_than_overrunning():
    """The cap wins over max_decoys; underfill is the visible consequence."""
    with tempfile.TemporaryDirectory() as tmp:
        result = _run(_build(tmp))
    total = sum(len(d) for _, d in result)
    assert total <= len(POOL_EXTRA) * CAP
    assert any(len(d) < MAX_DECOYS for _, d in result)
