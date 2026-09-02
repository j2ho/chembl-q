# chembl_curator/labeler.py

"""Label (target, compound) pairs as active, inactive, or unlabelled.

Pure logic: takes measurement records and returns labels. No database access,
so the rules can be unit tested independently of the ChEMBL schema.

Label rules, applied per (target, compound) pair over all its measurements:

    active              at least one exact potency measurement <= active_max_nm
    inactive_annotated  at least one explicit inactive activity_comment
    inactive_potency    at least one censored measurement (">") at or above
                        inactive_min_nm, i.e. tested that high with no binding
    weak                at least one exact potency measurement above
                        active_max_nm (binding observed, but too weak)
    uninformative       only censored measurements below inactive_min_nm
                        (tested too low to say anything)
    conflict            evidence for both active and inactive, close enough to
                        the active cutoff that neither side wins

When a pair carries both kinds of evidence, potency decides. A compound whose
best exact measurement sits a full log unit below the active cutoff is an
unambiguous binder and the stray inactive call is treated as noise: it stays
active, tagged "contested". Near the cutoff the disagreement is real and the
pair is dropped from both label sets.

Supporting record counts are recorded but deliberately not used to decide.
ChEMBL record counts track publication habits, not evidence quality; a single
paper can deposit dozens of replicates of one measurement.

Only "exact" relations carry binding evidence; a censored ">" record says a
compound was not seen to bind up to the tested concentration, which is
non-binding evidence whose strength is set by that concentration. Exact
measurements are therefore never treated as inactive evidence regardless of
magnitude.

Values are only interpreted when the measurement reports a potency type in
nM or uM. Annotated inactives are frequently reported as "% Control",
"Inhibition" or "Activity", whose numbers are not concentrations.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, FrozenSet, Iterable, List, Optional, Tuple

# Depositor activity calls that mean "this compound did not work".
# Compared case-insensitively. Deliberately excludes "inconclusive",
# "not determined" and free-text percentage remarks.
INACTIVE_COMMENTS = frozenset({"inactive", "not active", "no inhibition"})

# Potency measurement types whose standard_value is a concentration.
# Must stay in step with CurationConfig.activity_types: the active query and
# these rules have to agree on what counts as a potency measurement, or a
# compound can be missing from actives.tsv yet still be excluded from the
# inactive set. Callers should pass config.activity_types explicitly.
POTENCY_TYPES = frozenset({"Kd", "Ki", "IC50", "EC50"})

# Relations that assert a measured binding event.
EXACT_RELATIONS = frozenset({"=", "<=", "<"})
# Relations that assert no binding was seen up to the reported value.
CENSORED_RELATIONS = frozenset({">", ">=", ">>"})

UNIT_TO_NM = {"nM": 1.0, "uM": 1000.0}


class Label(str, Enum):
    ACTIVE = "active"
    INACTIVE_ANNOTATED = "inactive_annotated"
    INACTIVE_POTENCY = "inactive_potency"
    WEAK = "weak"
    UNINFORMATIVE = "uninformative"
    CONFLICT = "conflict"


@dataclass(frozen=True)
class Measurement:
    """One ChEMBL activity row, already mapped to UniProt and ChEMBL compound id."""

    target: str
    compound: str
    std_type: Optional[str]
    relation: Optional[str]
    value: Optional[float]
    units: Optional[str]
    pchembl: Optional[float]
    comment: Optional[str]


@dataclass
class PairVerdict:
    """Outcome for one (target, compound) pair."""

    label: Label
    evidence: Tuple[str, ...]
    best_pchembl: Optional[float]
    best_active_nm: Optional[float]
    n_active_records: int = 0
    n_inactive_records: int = 0

    @property
    def is_inactive(self) -> bool:
        return self.label in (Label.INACTIVE_ANNOTATED, Label.INACTIVE_POTENCY)

    @property
    def is_contested(self) -> bool:
        """Kept active despite an inactive call, on potency grounds."""
        return self.label is Label.ACTIVE and "contested" in self.evidence


def concentration_nm(
    m: Measurement,
    potency_types: FrozenSet[str] = POTENCY_TYPES,
) -> Optional[float]:
    """Concentration in nM, or None if this row does not report one."""
    if m.value is None or m.std_type not in potency_types:
        return None
    factor = UNIT_TO_NM.get(m.units or "")
    if factor is None:
        return None
    return float(m.value) * factor


def is_inactive_comment(comment: Optional[str]) -> bool:
    return bool(comment) and comment.strip().lower() in INACTIVE_COMMENTS


def classify_pair(
    measurements: Iterable[Measurement],
    active_max_nm: float = 10_000.0,
    inactive_min_nm: float = 100_000.0,
    conflict_decisive_nm: float = 1_000.0,
    potency_types: FrozenSet[str] = POTENCY_TYPES,
) -> PairVerdict:
    """Reduce every measurement of one (target, compound) pair to a single label.

    conflict_decisive_nm: an active measurement at or below this potency
    overrides a contradicting inactive call. Defaults to a tenth of the
    active cutoff.
    """
    has_active = has_annotated = has_censored_inactive = False
    has_weak = has_uninformative = False
    n_active = n_inactive = 0
    best_pchembl: Optional[float] = None
    best_active_nm: Optional[float] = None

    for m in measurements:
        if is_inactive_comment(m.comment):
            has_annotated = True
            n_inactive += 1

        conc = concentration_nm(m, potency_types)
        if conc is None:
            continue

        if m.relation in EXACT_RELATIONS:
            if conc <= active_max_nm:
                has_active = True
                n_active += 1
                if best_active_nm is None or conc < best_active_nm:
                    best_active_nm = conc
                if m.pchembl is not None:
                    p = float(m.pchembl)
                    if best_pchembl is None or p > best_pchembl:
                        best_pchembl = p
            else:
                has_weak = True
        elif m.relation in CENSORED_RELATIONS:
            if conc >= inactive_min_nm:
                has_censored_inactive = True
                n_inactive += 1
            else:
                has_uninformative = True

    evidence: List[str] = []
    if has_annotated:
        evidence.append("annotated")
    if has_censored_inactive:
        evidence.append("potency")

    def verdict(label: Label, ev: Tuple[str, ...], keep_potency: bool = True) -> PairVerdict:
        return PairVerdict(
            label, ev,
            best_pchembl if keep_potency else None,
            best_active_nm if keep_potency else None,
            n_active, n_inactive,
        )

    if has_active and (has_annotated or has_censored_inactive):
        decisive = best_active_nm is not None and best_active_nm <= conflict_decisive_nm
        if decisive:
            return verdict(Label.ACTIVE, ("contested",) + tuple(evidence))
        return verdict(Label.CONFLICT, tuple(evidence))
    if has_active:
        return verdict(Label.ACTIVE, ())
    if has_annotated:
        return verdict(Label.INACTIVE_ANNOTATED, tuple(evidence), keep_potency=False)
    if has_censored_inactive:
        return verdict(Label.INACTIVE_POTENCY, tuple(evidence), keep_potency=False)
    if has_weak:
        return verdict(Label.WEAK, (), keep_potency=False)
    return verdict(Label.UNINFORMATIVE, (), keep_potency=False)


def classify_all(
    measurements: Iterable[Measurement],
    active_max_nm: float = 10_000.0,
    inactive_min_nm: float = 100_000.0,
    conflict_decisive_nm: float = 1_000.0,
    potency_types: FrozenSet[str] = POTENCY_TYPES,
) -> Dict[Tuple[str, str], PairVerdict]:
    """Group measurements by (target, compound) and classify each group."""
    grouped: Dict[Tuple[str, str], List[Measurement]] = {}
    for m in measurements:
        grouped.setdefault((m.target, m.compound), []).append(m)
    return {
        key: classify_pair(rows, active_max_nm, inactive_min_nm,
                           conflict_decisive_nm, potency_types)
        for key, rows in grouped.items()
    }
