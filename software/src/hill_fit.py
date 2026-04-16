"""R10 Hill fit (reparametrized, weighted NLS) + R11 weighted R².

Hill model:
    y(x) = baseline + (top - baseline) * x^n / (K_D^n + x^n)

Reparametrization for numerical stability:
    amplitude = log(top - baseline)  =>  top = baseline + exp(amplitude)

scipy curve_fit uses sigma_j = 1/sqrt(w_j) (equivalent to weighted least-squares).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from constants import DELTA


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


def _hill_3p(x: np.ndarray, log_kd: float, n: float, amplitude: float, baseline_fixed: float) -> np.ndarray:
    """Hill evaluated with baseline pinned; returns y values."""
    kd = math.exp(log_kd)
    top_minus_base = math.exp(amplitude)
    xn = np.power(x, n)
    return baseline_fixed + top_minus_base * xn / (kd**n + xn)


def _hill_4p(x: np.ndarray, log_kd: float, n: float, amplitude: float, baseline: float) -> np.ndarray:
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
    y_min, y_max = float(np.min(y)), float(np.max(y))
    top_guess = max(y_max, baseline + 1e-6)
    amp_guess = math.log(max(top_guess - baseline, 1e-6))
    half = baseline + (top_guess - baseline) / 2.0
    # Concentration where signal first exceeds the midpoint — coarse K_D estimate.
    above = np.where(y >= half)[0]
    if above.size > 0 and x[above[0]] > 0:
        kd_guess = float(x[above[0]])
    else:
        kd_guess = float(np.median(x[x > 0])) if np.any(x > 0) else 1.0
    return math.log(max(kd_guess, 1e-12)), 1.0, amp_guess


def fit_one_clonotype(
    x: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    baseline_fixed: float | None,
    *,
    bin_mode: bool,
    max_bin_label: int | None,
    kd_lo: float = 1e-15,
    kd_hi: float = 1e3,
    n_lo: float = 0.1,
    n_hi: float = 10.0,
) -> FitResult:
    """Single-clonotype Hill fit. Pure numpy inputs; no polars involvement.

    Returns FitResult with `converged=False, reason="convergence_failure"` on:
      - scipy.curve_fit raising
      - y dynamic range < DELTA after fit (top − baseline < δ)
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)

    # sigma = 1/sqrt(w); guard against zero weights by treating them as very uncertain.
    sigma = np.where(w > 0, 1.0 / np.sqrt(w), 1e12)

    amp_hi_mode = (
        math.log(max_bin_label - 1) if bin_mode and max_bin_label and max_bin_label > 1
        else math.log(0.95)
    )
    # Keep the lower bound well below log(DELTA) so the solver can land below δ
    # for truly flat signals; the `top - baseline < DELTA` post-fit check rejects them.
    amp_lo = math.log(DELTA * 1e-6)

    try:
        if baseline_fixed is not None:
            log_kd0, n0, amp0 = _initial_guesses(x, y, baseline_fixed)
            # Clip initial guesses into bounds.
            log_kd0 = min(max(log_kd0, math.log(kd_lo)), math.log(kd_hi))
            n0 = min(max(n0, n_lo), n_hi)
            amp0 = min(max(amp0, amp_lo), amp_hi_mode)

            def model(x_, log_kd, n, amp):
                return _hill_3p(x_, log_kd, n, amp, baseline_fixed)

            popt, _ = curve_fit(
                model,
                x,
                y,
                p0=[log_kd0, n0, amp0],
                sigma=sigma,
                absolute_sigma=False,
                bounds=(
                    [math.log(kd_lo), n_lo, amp_lo],
                    [math.log(kd_hi), n_hi, amp_hi_mode],
                ),
                method="trf",
                maxfev=10_000,
            )
            log_kd, n, amp = popt
            baseline_out = baseline_fixed
            y_hat = _hill_3p(x, log_kd, n, amp, baseline_fixed)
        else:
            y_min = float(np.min(y))
            log_kd0, n0, amp0 = _initial_guesses(x, y, y_min)
            log_kd0 = min(max(log_kd0, math.log(kd_lo)), math.log(kd_hi))
            n0 = min(max(n0, n_lo), n_hi)
            amp0 = min(max(amp0, amp_lo), amp_hi_mode)
            base_lo, base_hi = 0.0, float(max_bin_label - 1) if bin_mode and max_bin_label else 1.0
            base0 = min(max(y_min, base_lo), base_hi)

            popt, _ = curve_fit(
                _hill_4p,
                x,
                y,
                p0=[log_kd0, n0, amp0, base0],
                sigma=sigma,
                absolute_sigma=False,
                bounds=(
                    [math.log(kd_lo), n_lo, amp_lo, base_lo],
                    [math.log(kd_hi), n_hi, amp_hi_mode, base_hi],
                ),
                method="trf",
                maxfev=10_000,
            )
            log_kd, n, amp, baseline_out = popt
            y_hat = _hill_4p(x, log_kd, n, amp, baseline_out)

    except (RuntimeError, ValueError, OptimizeWarning):
        return FitResult(
            kd=None, n=None, top=None, baseline=None,
            r2_w=None, y_hat=None, converged=False, reason="convergence_failure",
        )

    top_minus_base = math.exp(amp)
    if top_minus_base < DELTA:
        return FitResult(
            kd=None, n=None, top=None, baseline=None,
            r2_w=None, y_hat=None, converged=False, reason="convergence_failure",
        )

    kd = math.exp(log_kd)
    r2 = weighted_r2(y, y_hat, w)
    return FitResult(
        kd=kd,
        n=float(n),
        top=baseline_out + top_minus_base,
        baseline=float(baseline_out),
        r2_w=r2,
        y_hat=y_hat,
        converged=True,
        reason=None,
    )
