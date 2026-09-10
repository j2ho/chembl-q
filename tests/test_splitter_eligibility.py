#!/usr/bin/env python3
"""The test split may never contain a target that clusters with PDBbind/BioLiP."""

from pathlib import Path

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


def _splitter(demoted=()):
    s = TargetSplitter(valid_frac=1.0, log_level="CRITICAL")
    s.demoted_targets = set(demoted)
    return s


def test_demoted_target_is_labelled_train():
    """Cluster eligibility alone does not make a target test.

    P1 sits in a ChEMBL-only cluster, so valid_reps says test, but it was
    demoted for direct identity to a train sequence.
    """
    s = _splitter(demoted={"P1"})
    splits = s._assign_splits(
        chembl_actives={"P1": ["c1"], "P2": ["c2"]},
        member_to_rep={"chembl.P1": "r1", "chembl.P2": "r2"},
        valid_reps={"r1", "r2"},
    )
    assert splits["P1"] == "train"
    assert splits["P2"] == "test"


def test_without_demotion_the_same_target_is_test():
    """Guard for the test above: P1 must be test when nothing demotes it,
    otherwise the assertion would pass for the wrong reason."""
    s = _splitter()
    splits = s._assign_splits(
        chembl_actives={"P1": ["c1"]},
        member_to_rep={"chembl.P1": "r1"},
        valid_reps={"r1"},
    )
    assert splits["P1"] == "test"


def test_summary_cannot_disagree_with_the_split_files():
    """chembl_targets.tsv is written from the assignment, not from a second
    reading of valid_reps. Passing a label the cluster state contradicts must
    produce that label, which is what makes the two files agree by
    construction."""
    import inspect

    sig = inspect.signature(TargetSplitter._write_chembl_targets)
    assert "target_split" in sig.parameters, \
        "_write_chembl_targets must be handed the assignment"
    assert "valid_reps" not in sig.parameters, \
        "_write_chembl_targets must not be able to recompute the split"


def test_external_coverage_is_measured_over_the_external_sequence():
    """A crystallised domain covers little of a large UniProt entry.

    Requiring 80% of the ChEMBL query let 120 test targets through whose own
    representative PDB entry was sitting in BioLiP or PDBbind at pocket RMSD
    0.000; the worst covered 7% of the query (IGF2R, 2,491 aa against a 182 aa
    construct) while covering the whole external sequence.
    """
    import tempfile
    from unittest.mock import patch

    rows = [
        # query=external, target=chembl, fident, alnlen, qlen(ext), tlen(chembl)
        ("pdbbind.6n5x", "chembl.P11717", "0.99", "180", "182", "2491"),
        ("biolip.1abc_LIG_A_1", "chembl.P99999", "0.99", "60", "300", "100"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        (td / "external_hits.tsv").write_text(
            "".join("\t".join(r) + "\n" for r in rows))

        s = TargetSplitter(seqid=0.4, log_level="CRITICAL")
        with patch("chembl_curator.splitter.subprocess.run"):
            blocked = s._find_external_homologues(td / "q.fasta",
                                                  td / "e.fasta", td)

    assert "P11717" in blocked, \
        "a fully covered external domain inside a large target must block it"
    assert "P99999" not in blocked, \
        "a partial external sequence is a fragment hit and must not block"


def test_external_coverage_threshold_is_configurable():
    import tempfile
    from unittest.mock import patch

    # query=external (100 aa), target=chembl; 50 aligned is half the external
    row = ("biolip.1abc_LIG_A_1", "chembl.P1", "0.99", "50", "100", "1000")
    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        (td / "external_hits.tsv").write_text("\t".join(row) + "\n")

        strict = TargetSplitter(seqid=0.4, external_coverage=0.8,
                                log_level="CRITICAL")
        loose = TargetSplitter(seqid=0.4, external_coverage=0.5,
                               log_level="CRITICAL")
        with patch("chembl_curator.splitter.subprocess.run"):
            assert "P1" not in strict._find_external_homologues(
                td / "q.fasta", td / "e.fasta", td)
            assert "P1" in loose._find_external_homologues(
                td / "q.fasta", td / "e.fasta", td)


def _targets_fixture(tmp, with_external):
    """A data dir holding one train and one test target, each with decoys."""
    data = Path(tmp) / "data"
    for u in ("P1", "P2"):
        (data / u).mkdir(parents=True)
        (data / u / "decoys.tsv").write_text("active\tdecoys\nc1\td1;d2\n")
    if with_external:
        (data / "external_pocket_best.tsv").write_text(
            "chembl_target\texternal_entry\tpocket_rmsd\tn_matched"
            "\tn_pocket_chembl\tn_pocket_external\n"
            "P1\tbiolip.1abc_LIG_A_1\t0.372\t28\t53\t32\n")
    return data


def test_external_pocket_columns_are_reported_not_filtered():
    """Stage 7's verdict rides along in the summary; it must not move a target.

    Cutting on pocket similarity removes the data-rich targets, so the number
    is published and the reader decides.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        data = _targets_fixture(tmp, with_external=True)
        s = TargetSplitter(log_level="CRITICAL")
        out = s._write_chembl_targets(
            data, Path(tmp), {"P1": ["c1"], "P2": ["c1"]},
            {"P1": "test", "P2": "test"})

        # rstrip("\n") only: strip() would eat the trailing empty columns of
        # the last row, which are exactly what this test checks.
        rows = [l.split("\t") for l in out.read_text().rstrip("\n").split("\n")]
    header, p1, p2 = rows[0], rows[1], rows[2]
    assert header[-2:] == ["ext_pocket_rmsd", "ext_pocket_entry"]
    assert p1[1] == "test", "a close external pocket must not demote the target"
    assert p1[-2] == "0.372" and p1[-1] == "biolip.1abc_LIG_A_1"
    assert p2[-2] == "" and p2[-1] == "", \
        "a target with no external match must be blank, not zero"


def test_columns_are_absent_when_stage_8_has_not_run():
    """Blank would be honest, but a 0.0 default would read as 'identical to an
    external pocket'. Leaving the columns off says 'not measured'."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        data = _targets_fixture(tmp, with_external=False)
        s = TargetSplitter(log_level="CRITICAL")
        out = s._write_chembl_targets(
            data, Path(tmp), {"P1": ["c1"]}, {"P1": "train"})
        header = out.read_text().split("\n")[0].split("\t")
    assert header == ["uniprot", "split", "n_actives", "n_decoys"]


def test_owning_an_external_pdb_entry_blocks_a_target():
    """Sequence search cannot be relied on to catch this.

    IGF2R is 2,491 residues and its PDBbind entry 6N5X is a 182-residue
    construct; mmseqs reports a 31-residue alignment at 100% identity, which
    clears no coverage rule in either search direction. The PDB identifier
    settles it.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        ext = td / "external.fasta"
        ext.write_text(">pdbbind.6n5x\nAAAA\n>biolip.1abc_LIG_A_1\nCCCC\n")
        data = td / "data"
        for target, pdbs in (("P11717", "6n5x 1xyz"),       # owns 6N5X
                             ("P00001", "1ABC 2def"),        # owns 1abc, case differs
                             ("P00002", "9zzz 8yyy")):       # owns neither
            (data / target).mkdir(parents=True)
            (data / target / "pdbid.list").write_text(pdbs.replace(" ", "\n") + "\n")

        s = TargetSplitter(log_level="CRITICAL")
        blocked = s._find_external_structure_sharing(data, ext)

    assert blocked == {"P11717", "P00001"}, blocked


def test_structure_sharing_ignores_targets_with_no_pdb_list():
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        ext = td / "external.fasta"
        ext.write_text(">pdbbind.6n5x\nAAAA\n")
        data = td / "data"
        (data / "P00003").mkdir(parents=True)      # no pdbid.list at all
        s = TargetSplitter(log_level="CRITICAL")
        assert s._find_external_structure_sharing(data, ext) == set()


def test_homologue_search_reads_coverage_from_the_query_side():
    """The external sequence is the query, so its coverage is alnlen/qlen.

    Reading tlen instead would measure coverage of the ChEMBL target, which is
    the bug the external-coverage fix removed in the first place.
    """
    import tempfile
    from unittest.mock import patch

    rows = [
        # query=external, target=chembl, fident, alnlen, qlen(ext), tlen(chembl)
        ("pdbbind.6n5x", "chembl.P11717", "0.99", "180", "182", "2491"),
        ("biolip.1abc_LIG_A_1", "chembl.P99999", "0.99", "60", "300", "100"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        (td / "external_hits.tsv").write_text(
            "".join("\t".join(r) + "\n" for r in rows))
        s = TargetSplitter(seqid=0.4, log_level="CRITICAL")
        with patch("chembl_curator.splitter.subprocess.run"):
            blocked = s._find_external_homologues(td / "c.fasta", td / "e.fasta", td)

    assert "P11717" in blocked, "a fully covered external sequence must block"
    assert "P99999" not in blocked, "a fragment of the external sequence must not"


def test_demotion_measures_coverage_on_both_sides():
    """A short train sequence sitting whole inside a long test target counts.

    Q08209 is 521 residues and the train target P67775 is 309; their 283-residue
    alignment at 0.43 identity covers 92% of the train sequence and 54% of the
    test one. Measuring only over the query let it stay in the test set, which
    is the same coverage-direction mistake the external rule already fixed.
    """
    import tempfile
    from unittest.mock import patch

    rows = [
        # query=test target, target=train side, fident, alnlen, qlen, tlen
        ("Q08209", "chembltrain.P67775", "0.43", "283", "521", "309"),
        ("P99998", "chembltrain.P00001", "0.99", "40", "500", "600"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        train_fasta = td / "train.fasta"
        train_fasta.write_text(">chembltrain.P67775\nMMMM\n")
        # round 2 finds nothing, so the loop stops after one pass
        (td / "demote_hits1.tsv").write_text(
            "".join("\t".join(r) + "\n" for r in rows))
        (td / "demote_hits2.tsv").write_text("")

        s = TargetSplitter(seqid=0.4, log_level="CRITICAL")
        with patch("chembl_curator.splitter.subprocess.run"):
            demoted = s._demote_test_targets_close_to_train(
                {"Q08209", "P99998"}, train_fasta,
                {"Q08209": "AAAA", "P99998": "CCCC"}, td)

    assert "Q08209" in demoted, \
        "a train sequence covered end to end must demote the test target"
    assert "P99998" not in demoted, \
        "40 residues of either sequence is a fragment, not evidence"
