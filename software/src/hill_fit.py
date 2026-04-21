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
KD_LO: float = 1e-15
KD_HI: float = 1e3
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


def _failure() -> FitResult:
    return FitResult(
        kd=None,
        n=None,
        top=None,
        baseline=None,
        r2_w=None,
        y_hat=None,
        converged=False,
        reason="convergence_failure",
    )


def fit_one_clonotype(
    x: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    baseline_fixed: float | None,
    *,
    bin_mode: bool,
    max_bin_label: int | None,
) -> FitResult:
    """Single-clonotype Hill fit. Pure numpy inputs; no polars involvement.

    Returns FitResult with `converged=False, reason="convergence_failure"` on:
      - scipy.curve_fit raising
      - top − baseline < mode-specific δ after fit (bin: 0.5, no-bin: 0.05)
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)

    # sigma = 1/sqrt(w); guard against zero weights by treating them as very uncertain.
    sigma = np.where(w > 0, 1.0 / np.sqrt(w), 1e12)

    delta = DELTA_BIN if bin_mode else DELTA_NO_BIN
    amp_hi = _amplitude_upper_bound(bin_mode, max_bin_label)
    amp_lo = math.log(delta)

    baseline_known = baseline_fixed is not None
    anchor_for_guesses = baseline_fixed if baseline_known else float(np.min(y))
    log_kd0, n0, amp0 = _initial_guesses(x, y, anchor_for_guesses)
    log_kd0 = min(max(log_kd0, math.log(KD_LO)), math.log(KD_HI))
    n0 = min(max(n0, N_LO), N_HI)
    amp0 = min(max(amp0, amp_lo), amp_hi)

    if baseline_known:
        p0 = [log_kd0, n0, amp0]
        lo = [math.log(KD_LO), N_LO, amp_lo]
        hi = [math.log(KD_HI), N_HI, amp_hi]

        def model(x_, log_kd, n, amp):
            return _hill(x_, log_kd, n, amp, baseline_fixed)
    else:
        # R10 baseline bounds — signal-mode dependent, δ-aware so `top ≥ baseline + δ`
        # is always feasible within top's upper bound:
        #   bin:    baseline ∈ [1,   max_bin_label − 0.5]
        #   no-bin: baseline ∈ [0,   0.95]
        if bin_mode and max_bin_label:
            base_lo = 1.0
            base_hi = float(max_bin_label) - 0.5
        else:
            base_lo = 0.0
            base_hi = 0.95
        base0 = min(max(anchor_for_guesses, base_lo), base_hi)
        p0 = [log_kd0, n0, amp0, base0]
        lo = [math.log(KD_LO), N_LO, amp_lo, base_lo]
        hi = [math.log(KD_HI), N_HI, amp_hi, base_hi]
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

    if baseline_known:
        log_kd, n, amp = popt
        baseline_out = baseline_fixed
    else:
        log_kd, n, amp, baseline_out = popt
    y_hat = _hill(x, log_kd, n, amp, baseline_out)

    # A fit pinned at amp_lo (within FP tolerance) means the data wanted a
    # dynamic range below δ; spec says reject as convergence_failure.
    if amp <= amp_lo + 1e-6:
        return _failure()
    top_minus_base = math.exp(amp)

    return FitResult(
        kd=math.exp(log_kd),
        n=float(n),
        top=float(baseline_out) + top_minus_base,
        baseline=float(baseline_out),
        r2_w=weighted_r2(y, y_hat, w),
        y_hat=y_hat,
        converged=True,
        reason=None,
    )
