"""R10 Hill fit (reparametrized, weighted NLS) + R11 weighted R².

Hill model:
    y(x) = baseline + (top - baseline) * x^n / (Kd^n + x^n)

Reparametrization for numerical stability:
    amplitude = log(top - baseline)  =>  top = baseline + exp(amplitude)

scipy curve_fit uses sigma_j = 1/sqrt(w_j) (equivalent to weighted least-squares).
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from constants import DELTA_BIN, DELTA_NO_BIN

# Parameter bounds for the Hill fit — tuned once, not caller-configurable.
# Kd may sit this many decades beyond the titrated range, in the data's own concentration unit.
KD_RANGE_DECADES: float = 3.0
N_LO: float = 0.1
N_HI: float = 10.0


@dataclass
class FitResult:
    kd: float | None
    n: float | None
    top: float | None
    baseline: float | None
    r2_w: float | None
    y_hat: np.ndarray | None
    converged: bool
    reason: str | None  # "convergence_failure" on internal failures; None on success


def _hill(x: np.ndarray, log_kd: float, n: float, amplitude: float, baseline: float) -> np.ndarray:
    kd = math.exp(log_kd)
    top_minus_base = math.exp(amplitude)
    xn = np.power(x, n)
    return baseline + top_minus_base * xn / (kd**n + xn)


def weighted_r2(y: np.ndarray, y_hat: np.ndarray, w: np.ndarray) -> float:
    """R11: R²_w = 1 − Σ(w·(y−ŷ)²) / Σ(w·(y−ȳ_w)²), ȳ_w = Σ(w·y)/Σ(w).

    Zero-weight points contribute nothing to numerator or denominator. Negative
    values are NOT clamped (spec explicit). Returns nan if Σw = 0 or if Σ(w·(y−ȳ_w)²) = 0.
    """
    w = np.asarray(w, dtype=float)
    y = np.asarray(y, dtype=float)
    y_hat = np.asarray(y_hat, dtype=float)
    sum_w = float(w.sum())
    if sum_w == 0.0:
        return float("nan")
    y_bar_w = float((w * y).sum() / sum_w)
    num = float((w * (y - y_hat) ** 2).sum())
    den = float((w * (y - y_bar_w) ** 2).sum())
    if den == 0.0:
        return float("nan")
    return 1.0 - num / den


def _initial_guesses(x: np.ndarray, y: np.ndarray, baseline: float) -> tuple[float, float, float]:
    """Return (log_kd0, n0, amplitude0). Robust to monotonic + flat-ish inputs."""
    _y_min, y_max = float(np.min(y)), float(np.max(y))
    top_guess = max(y_max, baseline + 1e-6)
    amp_guess = math.log(max(top_guess - baseline, 1e-6))
    half = baseline + (top_guess - baseline) / 2.0
    # Concentration where signal first exceeds the midpoint — coarse Kd estimate.
    above = np.where(y >= half)[0]
    if above.size > 0 and x[above[0]] > 0:
        kd_guess = float(x[above[0]])
    else:
        kd_guess = float(np.median(x[x > 0])) if np.any(x > 0) else 1.0
    return math.log(max(kd_guess, 1e-12)), 1.0, amp_guess


def _amplitude_upper_bound(bin_mode: bool, max_bin_label: int | None) -> float:
    """R10: top ≤ max_bin_label - 1 in bin mode; ≤ 1 (freq) in no-bin mode."""
    if bin_mode and max_bin_label and max_bin_label > 1:
        return math.log(max_bin_label - 1)
    return math.log(0.95)


def _failure(reason: str = "convergence_failure") -> FitResult:
    return FitResult(
        kd=None,
        n=None,
        top=None,
        baseline=None,
        r2_w=None,
        y_hat=None,
        converged=False,
        reason=reason,
    )


def fit_one_clonotype(
    x: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    baseline_fixed: float | None,
    *,
    bin_mode: bool,
    max_bin_label: int | None,
    pin_n: bool = False,
    amp_prior_q: float | None = None,
) -> FitResult:
    """Single-clonotype Hill fit. Pure numpy inputs; no polars involvement.

    `pin_n` fixes the Hill exponent at 1, as the Adams 2016 model does. `amp_prior_q` raises the
    amplitude lower bound to that fraction of the bin range; None leaves it at δ.

    The amplitude floor bounds the optimiser. A fit that still lands on it had no dynamic range to
    measure, and is reported as `insufficient_dynamic_range` rather than as a bad fit.

    Returns FitResult with `converged=False` on:
      - scipy.curve_fit raising, or a flat signal with an undefined weighted R²
        (reason "convergence_failure")
      - a fit sitting on the amplitude floor (reason "insufficient_dynamic_range")
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)

    # sigma = 1/sqrt(w); guard against zero weights by treating them as very uncertain.
    sigma = np.where(w > 0, 1.0 / np.sqrt(w), 1e12)

    delta = DELTA_BIN if bin_mode else DELTA_NO_BIN
    amp_hi = _amplitude_upper_bound(bin_mode, max_bin_label)
    amp_lo = math.log(delta)
    if amp_prior_q is not None and bin_mode and max_bin_label and max_bin_label > 1:
        span = float(max_bin_label) - 1.0
        # Plateau fraction q of the bin range implies a saturated mean bin of 1 + q*span.
        thr = span * amp_prior_q if baseline_fixed is None else (1.0 + span * amp_prior_q) - float(baseline_fixed)
        amp_lo = math.log(min(max(thr, 1e-6), math.exp(amp_hi) * 0.999))

    baseline_known = baseline_fixed is not None
    anchor_for_guesses = baseline_fixed if baseline_known else float(np.min(y))
    log_kd0, n0, amp0 = _initial_guesses(x, y, anchor_for_guesses)
    positive = x[x > 0]
    log_kd_lo = math.log(float(positive.min())) - KD_RANGE_DECADES * math.log(10.0)
    log_kd_hi = math.log(float(positive.max())) + KD_RANGE_DECADES * math.log(10.0)
    log_kd0 = min(max(log_kd0, log_kd_lo), log_kd_hi)
    n0 = min(max(n0, N_LO), N_HI)
    amp0 = min(max(amp0, amp_lo), amp_hi)

    if baseline_known and pin_n:
        p0 = [log_kd0, amp0]
        lo = [log_kd_lo, amp_lo]
        hi = [log_kd_hi, amp_hi]

        def model(x_: np.ndarray, log_kd: float, amp: float) -> np.ndarray:
            return _hill(x_, log_kd, 1.0, amp, baseline_fixed)
    elif baseline_known:
        p0 = [log_kd0, n0, amp0]
        lo = [log_kd_lo, N_LO, amp_lo]
        hi = [log_kd_hi, N_HI, amp_hi]

        def model(x_: np.ndarray, log_kd: float, n: float, amp: float) -> np.ndarray:
            return _hill(x_, log_kd, n, amp, baseline_fixed)
    else:
        # R10 baseline bounds — signal-mode dependent, δ-aware so `top ≥ baseline + δ`
        # is always feasible within top's upper bound:
        #   bin:    baseline ∈ [1,   max_bin_label − 0.5]
        #   no-bin: baseline ∈ [0,   0.95]
        if bin_mode and max_bin_label:
            base_lo = 1.0
            # Keep `top ≤ max_bin_label` reachable once the amplitude floor rises above δ.
            base_hi = max(base_lo, min(float(max_bin_label) - 0.5, float(max_bin_label) - math.exp(amp_lo)))
        else:
            base_lo = 0.0
            base_hi = 0.95
        base0 = min(max(anchor_for_guesses, base_lo), base_hi)
        if pin_n:
            p0 = [log_kd0, amp0, base0]
            lo = [log_kd_lo, amp_lo, base_lo]
            hi = [log_kd_hi, amp_hi, base_hi]

            def model(x_: np.ndarray, log_kd: float, amp: float, base: float) -> np.ndarray:
                return _hill(x_, log_kd, 1.0, amp, base)
        else:
            p0 = [log_kd0, n0, amp0, base0]
            lo = [log_kd_lo, N_LO, amp_lo, base_lo]
            hi = [log_kd_hi, N_HI, amp_hi, base_hi]
            model = _hill

    try:
        # Promote OptimizeWarning to an exception. scipy emits it (not raises) when
        # covariance can't be estimated or when bounds pin the solution — signals that
        # match the `_failure()` contract even though popt may otherwise be returned.
        with warnings.catch_warnings():
            warnings.filterwarnings("error", category=OptimizeWarning)
            popt, _ = curve_fit(
                model,
                x,
                y,
                p0=p0,
                sigma=sigma,
                absolute_sigma=False,
                bounds=(lo, hi),
                method="trf",
                maxfev=10_000,
            )
    except (RuntimeError, ValueError, OptimizeWarning):
        return _failure()

    if baseline_known and pin_n:
        log_kd, amp = popt
        n, baseline_out = 1.0, baseline_fixed
    elif baseline_known:
        log_kd, n, amp = popt
        baseline_out = baseline_fixed
    elif pin_n:
        log_kd, amp, baseline_out = popt
        n = 1.0
    else:
        log_kd, n, amp, baseline_out = popt
    y_hat = _hill(x, log_kd, n, amp, baseline_out)

    # A fit sitting on the amplitude floor has no dynamic range to measure: the assay could not
    # separate this clonotype from background at any concentration. That is a statement about the
    # data, not a verdict on the fit, so it is reported as its own reason rather than as a bad fit.
    #
    # Measured on the three Adams 2016 replicates at the default rung: admitting these rows instead
    # roughly doubles the table (1,706 -> 3,566 on replicate 1) and costs rank agreement
    # 0.909 -> 0.779, 0.889 -> 0.803 and 0.946 -> 0.863, buying 4, 3 and 0 of the 250 tightest.
    if amp <= amp_lo + 1e-6:
        return _failure("insufficient_dynamic_range")

    top_minus_base = math.exp(amp)
    r2 = weighted_r2(y, y_hat, w)
    # A flat signal has no variance to explain: R² is 0/0 and K_D is unconstrained.
    if math.isnan(r2):
        return _failure()

    return FitResult(
        kd=math.exp(log_kd),
        n=float(n),
        top=float(baseline_out) + top_minus_base,
        baseline=float(baseline_out),
        r2_w=r2,
        y_hat=y_hat,
        converged=True,
        reason=None,
    )
