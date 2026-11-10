"""R12: classification — affinity class + failure reason, one derivation per equation.

`classify` is the "hill_gate_index" derivation: v2.4.0's, plus the guard on an undefined R².
`classify_log_fluorescence` is the "log_fluorescence" derivation; see its docstring.
`rate_mean` rates a mean over synonymous variants or over replicates; see its docstring.

Spec truth table for `classify` (applied in order per row; boundaries inclusive on both sides of
the n-in-range check, inclusive on R² thresholds):

  Did not converge                      -> Failed, convergence_failure
  R² not finite (no signal variance)    -> Failed, low_r2
  R² < r2_threshold_failed              -> Failed, low_r2
  R² >= r2_threshold_good  and n in [nMin,nMax]    -> Good
  R² >= r2_threshold_good  and n out of range      -> Partial, n_out_of_range (downgrade)
  r2_threshold_failed <= R² < r2_threshold_good
     and n in [nMin,nMax]                -> Partial
     and n out of range                  -> Failed, n_out_of_range
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

from constants import KD_SD_GOOD, AffinityClass, FailureReason, FitParams, KdRangeSide

# A censored value's side of the titration: "below" the lowest concentration or "above" the highest.
Side = Literal["below", "above"]

# The value's place against the titrated range, as the kdRangeSide column reports it.
RANGE_SIDE: dict[Side | None, KdRangeSide] = {"below": "Below range", None: "In range", "above": "Above range"}

# Members censored at opposite walls with no side holding the majority.
REASON_CONFLICTING = "conflicting_censoring"


def censored_side(sides: list[Side]) -> Side | None:
    """The wall a value with no measured member sits at, or None when its members conflict.

    One side must hold at least half of the members and more than the other side. Averaging
    opposite walls would land inside the titrated range and mean nothing. The block's own rule: no
    paper defines one.
    """
    n = len(sides)
    count = {side: sides.count(side) for side in ("below", "above")}
    for side, other in (("below", "above"), ("above", "below")):
        if 2 * count[side] >= n and count[side] > count[other]:
            return side
    return None


def r2_undefined(r2: float | None) -> bool:
    """A signal with no variance gives weighted_r2 a zero denominator, so it returns NaN.

    The fitted curve is then an artefact of the amplitude bound, not a measurement: on a flat mean
    bin it reports a K_D of 28 M against a titration that stops at 1e-6. Such a row is Failed and
    carries no fitted numbers.
    """
    return r2 is not None and not math.isfinite(r2)


@dataclass
class Classification:
    affinity_class: AffinityClass
    failure_reason: FailureReason | None


def classify(
    r2: float | None,
    n: float | None,
    converged: bool,
    params: FitParams,
    *,
    precomputed_reason: FailureReason | None = None,
) -> Classification:
    """Return (class, reason). `precomputed_reason` is honoured when set (upstream gate)."""
    if precomputed_reason is not None:
        return Classification("Failed", precomputed_reason)
    if not converged or r2 is None or n is None:
        return Classification("Failed", "convergence_failure")

    n_in_range = params.n_min <= n <= params.n_max

    # Every comparison below reads False against NaN, so without this guard an undefined R² reaches
    # the n-in-range branch and reports Partial.
    if r2_undefined(r2):
        return Classification("Failed", "low_r2")
    if r2 < params.r2_threshold_failed:
        return Classification("Failed", "low_r2")
    if r2 >= params.r2_threshold_good:
        if n_in_range:
            return Classification("Good", None)
        return Classification("Partial", None)
    if n_in_range:
        return Classification("Partial", None)
    return Classification("Failed", "n_out_of_range")


def classify_log_fluorescence(
    *,
    range_side: KdRangeSide,
    log10_kd_se: float | None,
    r2: float | None,
    r2_gate_min: float | None,
    kd_se_gate_max: float | None = None,
) -> Classification:
    """The "log_fluorescence" derivation for a fit that converged. In order:

      Below range  K_D is pinned at the tight wall: it is at most the lowest concentration
      Above range  K_D is pinned at the weak wall: it is at least the highest concentration. A
                   curve that stayed flat and dim lands here, and so does one that rose without
                   reaching half-maximum
      Failed       low_r2: the R² gate is on and an in-range fit's R² is below it or undefined
      Failed       high_kd_se: an SE limit is set and an in-range fit's SE(log10 K_D) is above it
                   or undefined
      Partial      SE(log10 K_D) above KD_SD_GOOD or not finite
      Good         K_D inside the titrated range, SE(log10 K_D) at or below KD_SD_GOOD

    A bound is its own class rather than Partial, because which wall it sits at is what a reader
    acts on: below the range is the tightest binder in the library, above it is the weakest.

    Both gates read in-range fits only, as Phillips 2021 filters them (R² >= 0.8, SE <= 1). A bound has no curve
    through its midpoint for R² to measure. The Hill exponent gates nothing, so `n_out_of_range`
    never fires. The pipeline rates a gated clonotype or a fit that did not converge Failed
    before this runs.
    """
    if range_side != "In range":
        return Classification(range_side, None)
    if r2_gate_min is not None and (r2 is None or not math.isfinite(r2) or r2 < r2_gate_min):
        return Classification("Failed", "low_r2")
    if kd_se_gate_max is not None and (
        log10_kd_se is None or not math.isfinite(log10_kd_se) or log10_kd_se > kd_se_gate_max
    ):
        return Classification("Failed", "high_kd_se")
    if log10_kd_se is None or not math.isfinite(log10_kd_se) or log10_kd_se > KD_SD_GOOD:
        return Classification("Partial", None)
    return Classification("Good", None)


def bound_class(side: Side, own_classes: list[str], *, wall_classes: bool) -> str:
    """The class of a value censored at `side`, combined from fits whose own classes are `own_classes`.

    Under "log_fluorescence" a bound is rated by its wall, as a single fit is. "hill_gate_index"
    declares no wall class: a bound keeps the class its fits share, or is Partial when they differ.
    """
    if wall_classes:
        return RANGE_SIDE[side]
    return own_classes[0] if len(set(own_classes)) == 1 else "Partial"


# ---- A mean over members: the synonymous variants of one amino-acid variant, or the replicates of one
# variant. Both combinations share the value, the SE and the class rule below.

MemberState = Literal["measured", "below", "above", "failed"]
STATES: tuple[MemberState, ...] = ("measured", "below", "above", "failed")


@dataclass(frozen=True)
class Member:
    """One fit as it enters a mean: its state, its value at the wall or its own, its fit SE."""

    state: MemberState
    log10_kd: float | None
    log10_kd_se: float | None
    affinity_class: str


def member_from_fit(row: dict, log10_lo: float | None, log10_hi: float | None) -> Member:
    """A per-clonotype fit row as it enters a mean.

    A K_D flagged out of range enters at the nearer wall, with no SE. "log_fluorescence" pins its K_D
    at the wall already; "hill_gate_index" reports the raw value, so the midpoint of the range decides
    which wall it passed. A non-finite SE reads as none.
    """
    kd = row["kd"]
    cls = row["affinityClass"]
    if cls == "Failed" or kd is None or not math.isfinite(kd) or kd <= 0:
        return Member("failed", None, None, cls)
    value = math.log10(kd)
    if row["kdOutOfRange"] and log10_lo is not None and log10_hi is not None:
        if value <= (log10_lo + log10_hi) / 2:
            return Member("below", log10_lo, None, cls)
        return Member("above", log10_hi, None, cls)
    se = row["kdLog10Se"]
    return Member("measured", value, se if se is not None and math.isfinite(se) else None, cls)


def fit_error_term(entered: Sequence[Member]) -> float | None:
    """RMS of the fit SEs the members report, over sqrt(n); n counts every member in the mean.

    A member with no fit SE, a bound included, counts at the typical error of the others. None when
    no member reports a fit SE.
    """
    ses = [m.log10_kd_se for m in entered if m.log10_kd_se is not None]
    if not ses:
        return None
    return math.sqrt(math.fsum(se * se for se in ses) / len(ses)) / math.sqrt(len(entered))


def _majority(flags: list[bool]) -> bool:
    return 2 * sum(flags) > len(flags)


def rate_mean(
    entered: Sequence[Member], *, fits_report_se: bool
) -> tuple[float, float | None, float | None, Classification]:
    """The mean of `entered`, at least one of them measured: (log10 K_D, SE, SD, class).

    `fits_report_se` is True under the equation whose fits report an SE ("log_fluorescence"). It is
    read from the equation, not from the members: under "hill_gate_index" a synonymous combination
    gives its value a spread SE, and that value enters the mean over replicates.

    n = 1   the member's own SE; with none, its own class
    n >= 2  SE = max(sample SD / sqrt(n), fit_error_term), the spread alone when no member carries
            an SE; SD the sample SD, never floored
    Good    SE finite and at or below KD_SD_GOOD, and more than half the members measured; under an
            equation whose fits report no SE, or when no member carries an SE, also more than half
            of them rated Good. A member at a wall enters at the wall, so two of them agree exactly
            however wide the truth is
    Partial otherwise
    """
    n = len(entered)
    values = [m.log10_kd for m in entered]
    mean = math.fsum(values) / n
    no_fit_se = all(m.log10_kd_se is None for m in entered)
    if n == 1:
        se, sd = entered[0].log10_kd_se, None
        if se is None:
            return mean, None, None, Classification(entered[0].affinity_class, None)
    else:
        spread = statistics.stdev(values) / math.sqrt(n)
        fit_term = fit_error_term(entered)
        se = spread if fit_term is None else max(spread, fit_term)
        sd = statistics.stdev(values)
    good = math.isfinite(se) and se <= KD_SD_GOOD and _majority([m.state == "measured" for m in entered])
    # Reached at n = 1 only with an SE: under "hill_gate_index" a spread, so the member's own class rules.
    if not fits_report_se or no_fit_se:
        good = good and _majority([m.affinity_class == "Good" for m in entered])
    return mean, se, sd, Classification("Good" if good else "Partial", None)


@dataclass(frozen=True)
class Combined:
    """A mean over members, or the wall they share, and its class."""

    log10_kd: float | None
    log10_kd_se: float | None
    log10_kd_sd: float | None
    classification: Classification
    out_of_range: bool | None
    counts: dict[MemberState, int] = field(default_factory=dict)
    # The wall a bound sits at; None for a mean of measured values or for no value.
    side: Side | None = None

    @property
    def range_side(self) -> str | None:
        return None if self.log10_kd is None else RANGE_SIDE[self.side]


def failed(reason: FailureReason, counts: dict[MemberState, int] | None = None) -> Combined:
    """A combined value with no K_D, rated Failed for `reason`."""
    return Combined(None, None, None, Classification("Failed", reason), None, counts or {})


def combine_members(
    members: Sequence[Member],
    *,
    all_failed_reason: FailureReason,
    wall_classes: bool = True,
    fits_report_se: bool = True,
) -> Combined:
    """One value over `members`, and its class.

    Failed members are left out and counted. With no member left the value fails for
    `all_failed_reason`. With none of them measured the value is the wall they sit at
    (`censored_side`), or fails where they conflict. Otherwise it is their mean (`rate_mean`).
    `wall_classes` is True under the equation that rates a bound by its wall ("log_fluorescence");
    `fits_report_se` under the equation whose fits report an SE.
    """
    counts: dict[MemberState, int] = dict.fromkeys(STATES, 0)
    for m in members:
        counts[m.state] += 1
    entered = [m for m in members if m.state != "failed"]
    if not entered:
        return failed(all_failed_reason, counts)

    if not any(m.state == "measured" for m in entered):
        side = censored_side([m.state for m in entered])
        if side is None:
            return failed(REASON_CONFLICTING, counts)
        on_side = [m for m in entered if m.state == side]
        cls = bound_class(side, [m.affinity_class for m in on_side], wall_classes=wall_classes)
        return Combined(on_side[0].log10_kd, None, None, Classification(cls, None), True, counts, side)

    mean, se, sd, classification = rate_mean(entered, fits_report_se=fits_report_se)
    return Combined(mean, se, sd, classification, False, counts)
