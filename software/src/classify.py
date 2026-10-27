"""R12: classification — map (R²_w, n, converged) to affinity class + failure reason.

A fit that converged is always reported. R² and the Hill-coefficient range set the class; neither
discards the fit. An adequacy gate here cost 215/209/217 of the 250 tightest binders on the three
Adams 2016 replicates, measured by leave-one-out at the recommended configuration.

Truth table:

  Did not converge                                  -> Failed, convergence_failure
  R² >= r2_threshold_good and n in [nMin, nMax]     -> Good
  anything else that converged                      -> Partial
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from constants import AffinityClass, FailureReason, FitParams


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
    # A non-finite R² (a flat signal has no variance to explain) clears no threshold.
    if math.isfinite(r2) and r2 >= params.r2_threshold_good and n_in_range:
        return Classification("Good", None)
    return Classification("Partial", None)
