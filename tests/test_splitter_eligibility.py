#!/usr/bin/env python3
"""The test split may never contain a target that clusters with PDBbind/BioLiP."""

from chembl_curator.splitter import TargetSplitter


def _split(clusters, valid_frac=1.0, blocked=None):
    """valid_frac 1.0 is the shipped setting: every eligible cluster is test.

    Smaller fractions make the greedy balancer tie-break toward train on tiny
    inputs, which says nothing about eligibility.
    """
    s = TargetSplitter(valid_frac=valid_frac, log_level="CRITICAL")
    s.blocked_targets = set(blocked or ())
    return s._greedy_split(clusters)


def test_direct_homologue_blocks_even_a_pure_chembl_cluster():
    """Eligibility is per target, not per cluster."""
    clusters = {"r1": {"chembl.P1"}, "r2": {"chembl.P2"}}
    train, test = _split(clusters, blocked={"P1"})
    assert "r1" in train
    assert "r2" in test


def test_cluster_mate_without_a_homologue_is_not_punished():
    """P2 shares a cluster with a blocked target but has no external hit itself.

    Cluster-level blocking rejected 311 of 1,317 real targets this way.
    """
    clusters = {"r1": {"chembl.P1"}, "r2": {"chembl.P2"}, "r3": {"chembl.P3"}}
    train, test = _split(clusters, blocked={"P1"})
    assert {"r2", "r3"} <= test


def test_mixed_clusters_are_forced_to_train():
    clusters = {
        "r1": {"chembl.P1", "pdbbind.1abc"},
        "r2": {"chembl.P2", "biolip.2xyz_LIG_A_1"},
        "r3": {"chembl.P3"},
        "r4": {"chembl.P4"},
    }
    train, test = _split(clusters)
    assert "r1" in train and "r2" in train
    assert "r1" not in test and "r2" not in test
    assert test <= {"r3", "r4"}


def test_pure_chembl_clusters_can_reach_test():
    clusters = {f"r{i}": {f"chembl.P{i}"} for i in range(10)}
    train, test = _split(clusters)
    assert test, "eligible clusters must be able to reach test"
    assert not train & test


def test_blocked_chembl_targets_still_reach_train(capsys=None):
    """Forced-to-train clusters must be counted, not dropped from the tally.

    The tally used to run over the eligible subset only, so a run where every
    eligible cluster went to test logged "Actual train: chembl=0" while 539
    ChEMBL targets were in fact written to train.txt.
    """
    import io
    import logging

    stream = io.StringIO()
    handler = logging.StreamHandler(stream)
    logger = logging.getLogger("chembl_curator.splitter")
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    try:
        clusters = {
            "r1": {"chembl.P1", "pdbbind.1abc"},
            "r2": {"chembl.P2"},
        }
        s = TargetSplitter(valid_frac=1.0, log_level="INFO")
        s.blocked_targets = set()
        train, test = s._greedy_split(clusters)
    finally:
        logger.removeHandler(handler)

    assert "r1" in train and "r2" in test
    out = stream.getvalue()
    assert "Actual train: chembl=1" in out, out
    assert "chembl=0" not in out.split("Actual train:")[1].split("\n")[0]


def test_every_cluster_is_placed_exactly_once():
    clusters = {
        "r1": {"chembl.P1", "pdbbind.1abc"},
        "r2": {"chembl.P2"},
        "r3": {"chembl.P3"},
        "r4": {"biolip.9zzz_LIG_A_1"},
    }
    train, test = _split(clusters)
    assert train | test == set(clusters)
    assert not train & test


def test_external_only_cluster_never_lands_in_test():
    """A cluster with no ChEMBL member is not eligible either."""
    clusters = {
        "r1": {"pdbbind.1abc", "biolip.2xyz_LIG_A_1"},
        "r2": {"chembl.P2"},
    }
    train, test = _split(clusters)
    assert "r1" in train


def test_no_eligible_cluster_is_an_error_not_an_empty_test_set():
    clusters = {"r1": {"chembl.P1", "pdbbind.1abc"}}
    try:
        _split(clusters)
    except RuntimeError as e:
        assert "no test set" in str(e).lower() or "No ChEMBL-only" in str(e)
    else:
        raise AssertionError("silently produced an empty test set")


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
