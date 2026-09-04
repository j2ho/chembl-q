#!/usr/bin/env python3
"""Config flags to SQL conditions.

This used to re-type the curator's query-building by hand and print activity
counts against chembl_35.db, which meant it asserted nothing, silently did
nothing when the database was absent, and could drift from the code it was
supposed to be checking. It now exercises the real builder and needs no
database.
"""

from chembl_curator.config import CurationConfig
from chembl_curator.curator import ChEMBLCurator


def conditions(**kwargs):
    curator = ChEMBLCurator(config=CurationConfig(**kwargs), log_level="ERROR")
    return curator._shared_quality_conditions()


def test_flags_off_produce_no_conditions():
    assert conditions(require_standard_flag=False,
                      exclude_invalid_data=False,
                      exclude_duplicates=False) == []


def test_each_flag_contributes_its_own_condition():
    assert "a.standard_flag = 1" in conditions(require_standard_flag=True)
    assert "a.data_validity_comment IS NULL" in conditions(exclude_invalid_data=True)
    assert any("potential_duplicate" in c
               for c in conditions(exclude_duplicates=True))


def test_assay_quality_filters_are_rendered():
    got = conditions(min_confidence_score=8,
                     assay_types=["B"],
                     bao_formats=["BAO_0000357"])
    assert "ass.confidence_score >= 8" in got
    assert "ass.assay_type IN ('B')" in got
    assert "ass.bao_format IN ('BAO_0000357')" in got


def test_pchembl_floor_is_not_shared():
    """It is an active-only potency rule. A compound with no measurable
    affinity has no pChEMBL, so applying the floor to negatives would drop
    every one of them."""
    got = conditions(min_pchembl_value=5.0)
    assert not any("pchembl" in c.lower() for c in got)


def test_standard_type_is_not_shared():
    """The one filter that legitimately differs between the two queries: an
    inactive compound has no IC50 to report, so its result is filed under
    "% Control" or "Inhibition" instead."""
    got = conditions(activity_types=["Ki", "IC50"])
    assert not any("standard_type" in c for c in got)


def test_actives_and_negatives_share_the_same_base():
    """The point of the shared builder: neither query may quietly gain a
    quality filter the other lacks."""
    config = CurationConfig(require_standard_flag=True,
                            min_confidence_score=8,
                            bao_formats=["BAO_0000357"],
                            min_pchembl_value=5.0)
    curator = ChEMBLCurator(config=config, log_level="ERROR")
    base = curator._shared_quality_conditions()
    assert curator._shared_quality_conditions() == base
    assert len(base) == len(set(base)), "a condition is emitted twice"


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
    print("ok")
