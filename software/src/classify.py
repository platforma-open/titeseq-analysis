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

from constants import KD_SD_GOOD, AffinityClass, FailureReason, FitParams


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
    non_binder: bool,
    censored: bool,
    log10_kd_se: float | None,
    r2: float | None,
    r2_gate_min: float | None,
) -> Classification:
    """The "log_fluorescence" derivation for a fit that converged. In order:

      Non-binder  the curve is flat and dim, so its K_D is pinned at the weak wall
      Partial     the reported K_D is a bound: a flat bright curve, or a fit pinned at either wall
      Failed      low_r2: the R² gate is on and an in-range fit's R² is below it or undefined
      Partial     SE(log10 K_D) above KD_SD_GOOD or not finite
      Good        K_D inside the titrated range, SE(log10 K_D) at or below KD_SD_GOOD

    The R² gate reads in-range fits only, as Phillips 2021 filters them. A bound has no curve
    through its midpoint for R² to measure. The Hill exponent gates nothing, so `n_out_of_range`
    never fires. The pipeline rates a gated clonotype or a fit that did not converge Failed
    before this runs.
    """
    if non_binder:
        return Classification("Non-binder", None)
    if censored:
        return Classification("Partial", None)
    if r2_gate_min is not None and (r2 is None or not math.isfinite(r2) or r2 < r2_gate_min):
        return Classification("Failed", "low_r2")
    if log10_kd_se is None or not math.isfinite(log10_kd_se) or log10_kd_se > KD_SD_GOOD:
        return Classification("Partial", None)
    return Classification("Good", None)
