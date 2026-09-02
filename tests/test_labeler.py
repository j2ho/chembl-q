#!/usr/bin/env python3
"""Unit tests for the active/inactive labelling rules."""

from chembl_curator.labeler import (
    Label,
    Measurement,
    classify_all,
    classify_pair,
    concentration_nm,
)


def m(std_type="IC50", relation="=", value=None, units="nM", pchembl=None, comment=None):
    return Measurement("P00001", "CHEMBL1", std_type, relation, value, units,
                       pchembl, comment)


def test_exact_below_cutoff_is_active():
    v = classify_pair([m(value=500.0, pchembl=6.3)])
    assert v.label is Label.ACTIVE
    assert v.best_active_nm == 500.0
    assert v.best_pchembl == 6.3


def test_uM_units_are_converted():
    assert concentration_nm(m(value=5.0, units="uM")) == 5000.0
    assert classify_pair([m(value=5.0, units="uM")]).label is Label.ACTIVE
    assert classify_pair([m(value=50.0, units="uM")]).label is Label.WEAK


def test_exact_above_cutoff_is_weak_not_inactive():
    """An exact measurement is a binding event at any magnitude."""
    for value in (50_000.0, 150_000.0, 5_000_000.0):
        assert classify_pair([m(value=value)]).label is Label.WEAK


def test_censored_high_is_inactive():
    v = classify_pair([m(relation=">", value=100_000.0)])
    assert v.label is Label.INACTIVE_POTENCY
    assert v.evidence == ("potency",)


def test_censored_low_is_uninformative():
    """'>1 uM' says nothing about 10 uM, let alone 100 uM."""
    assert classify_pair([m(relation=">", value=1000.0)]).label is Label.UNINFORMATIVE
    assert classify_pair([m(relation=">", value=50_000.0)]).label is Label.UNINFORMATIVE


def test_annotated_inactive_without_value():
    v = classify_pair([m(std_type="% Control", value=None, comment="Not Active")])
    assert v.label is Label.INACTIVE_ANNOTATED
    assert v.evidence == ("annotated",)


def test_annotated_comment_matching_is_case_insensitive():
    for c in ("inactive", "Inactive", "NOT ACTIVE", "  no inhibition  "):
        assert classify_pair([m(comment=c)]).label is Label.INACTIVE_ANNOTATED


def test_inconclusive_is_not_inactive():
    for c in ("Inconclusive", "Not Determined", "ND", "Active"):
        assert classify_pair([m(comment=c)]).label is not Label.INACTIVE_ANNOTATED


def test_non_potency_values_are_not_concentrations():
    """'% Control = 5' must never be read as 5 nM."""
    assert concentration_nm(m(std_type="% Control", value=5.0)) is None
    assert classify_pair([m(std_type="Inhibition", value=5.0)]).label is Label.UNINFORMATIVE


def test_potency_types_follow_the_caller():
    """Must track config.activity_types, or the active query and these rules diverge."""
    only_ki = frozenset({"Ki"})
    rows = [m(std_type="IC50", value=500.0), m(comment="inactive")]
    assert classify_pair(rows).label is Label.ACTIVE
    # With IC50 out of scope the potency evidence disappears, so the pair is
    # inactive rather than silently falling out of both sets.
    assert classify_pair(rows, potency_types=only_ki).label is Label.INACTIVE_ANNOTATED


def test_borderline_active_plus_annotated_inactive_is_conflict():
    """Near the cutoff the disagreement is real, so drop the pair."""
    v = classify_pair([m(value=8600.0, pchembl=5.1), m(value=None, comment="inactive")])
    assert v.label is Label.CONFLICT
    assert "annotated" in v.evidence
    assert not v.is_inactive


def test_decisively_potent_beats_a_stray_inactive_call():
    """CHEMBL20 on P00918: 0.34 nM measured, one stray 'inactive' annotation."""
    v = classify_pair([m(value=0.34, pchembl=9.5), m(value=None, comment="inactive")])
    assert v.label is Label.ACTIVE
    assert v.is_contested
    assert v.best_active_nm == 0.34


def test_conflict_boundary_is_the_decisive_threshold():
    just_under = classify_pair([m(value=1000.0), m(comment="inactive")])
    just_over = classify_pair([m(value=1001.0), m(comment="inactive")])
    assert just_under.label is Label.ACTIVE
    assert just_over.label is Label.CONFLICT


def test_active_plus_censored_inactive_is_conflict():
    v = classify_pair([m(value=8600.0), m(relation=">", value=100_000.0)])
    assert v.label is Label.CONFLICT


def test_record_counts_are_reported_but_do_not_decide():
    """75 potent records vs 1 inactive call: potency decides, counts only inform."""
    rows = [m(value=8600.0) for _ in range(75)] + [m(comment="inactive")]
    v = classify_pair(rows)
    assert v.n_active_records == 75
    assert v.n_inactive_records == 1
    assert v.label is Label.CONFLICT


def test_annotated_plus_weak_is_not_conflict():
    """Inactive call alongside a 50 uM binding measurement is consistent."""
    v = classify_pair([m(value=50_000.0), m(value=None, comment="inactive")])
    assert v.label is Label.INACTIVE_ANNOTATED


def test_both_inactive_evidences_are_recorded():
    v = classify_pair([m(comment="inactive"), m(relation=">", value=100_000.0)])
    assert v.label is Label.INACTIVE_ANNOTATED
    assert v.evidence == ("annotated", "potency")


def test_classify_all_groups_by_pair():
    rows = [
        Measurement("P1", "C1", "IC50", "=", 100.0, "nM", 7.0, None),
        Measurement("P1", "C2", "IC50", ">", 100_000.0, "nM", None, None),
        Measurement("P2", "C1", "%", "=", 3.0, "%", None, "Not Active"),
    ]
    out = classify_all(rows)
    assert out[("P1", "C1")].label is Label.ACTIVE
    assert out[("P1", "C2")].label is Label.INACTIVE_POTENCY
    assert out[("P2", "C1")].label is Label.INACTIVE_ANNOTATED


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
