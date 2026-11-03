"""R12: classification — affinity class + failure reason, one derivation per equation.

`classify` is the "hill_gate_index" derivation: v2.4.0's, plus the guard on an undefined R².
`classify_log_fluorescence` is the "log_fluorescence" derivation; see its docstring.

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
from dataclasses import dataclass

from constants import KD_SD_GOOD, AffinityClass, FailureReason, FitParams, KdRangeSide


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
