"""R12: classification — affinity class + failure reason, one derivation per equation.

`classify` is the "hill_gate_index" derivation and is v2.4.0's, unchanged.
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

    # A signal with no variance gives a zero denominator, so weighted_r2 returns NaN. Every
    # comparison below reads False against NaN. Without this guard the row reaches the
    # n-in-range branch and reports Partial, carrying whatever K_D the fitter returned.
    if not math.isfinite(r2):
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
    converged: bool,
    censored: bool,
    log10_kd_se: float | None,
    precomputed_reason: FailureReason | None = None,
) -> Classification:
    """The "log_fluorescence" derivation.

      Failed   the fit did not converge, or too few reads or points to fit
      Partial  converged, but the reported K_D is a bound (censored, or pinned at a titration
               boundary), or SE(log10 K_D) is above KD_SD_GOOD or not finite
      Good     converged, K_D inside the titrated range, and SE(log10 K_D) at or below KD_SD_GOOD

    `precomputed_reason` is honoured when set (upstream gate). R² and the Hill exponent gate
    nothing here, so `low_r2` and `n_out_of_range` never fire.
    """
    if precomputed_reason is not None:
        return Classification("Failed", precomputed_reason)
    if not converged:
        return Classification("Failed", "convergence_failure")
    if censored:
        return Classification("Partial", None)
    if log10_kd_se is None or not math.isfinite(log10_kd_se) or log10_kd_se > KD_SD_GOOD:
        return Classification("Partial", None)
    return Classification("Good", None)
