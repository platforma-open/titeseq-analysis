"""The "log_fluorescence" equation: a Hill fit on mean log10 gate fluorescence.

Model, per clonotype:

    F̄(c) = log10( A * c / (c + K) + B )

fitted to the observed mean log10 fluorescence by weighted least squares, Hill exponent fixed at 1.

Gate numbers are proportional to mean *log* fluorescence, not to mean fluorescence, so fitting a
linear Hill curve to the gate index biases the result (Moulana 2022); measured on CR6261 that bias
is about 1.2 decades. Fitting in log space removes it.

Parameters are carried as log10 of the saturated top and of the baseline, so the amplitude
A = 10^top − 10^base is positive by construction and the plateau floor is a plain box bound on
`top`.

The module also carries the equation's own bound machinery — the gate scale, the observed-curve
class and the pin — because the rating reads whether the reported K_D is a bound.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from constants import GATE_RANGE_DECADES, GATE_WIDTH_DECADES, LF_RISE_THRESHOLD_DECADES
from fit_stats import weighted_r2

# K may sit this many decades beyond the titrated range before the optimiser is stopped. The value
# is a numerical guard, not a scientific claim: anything outside the range is reported as a bound.
KD_RANGE_DECADES: float = 3.0
# Smallest gap, in decades, kept between the baseline's upper bound and the plateau floor so the two
# box bounds can never cross.
_BOUND_GAP: float = 1e-3

# How the reported K_D relates to the true one when a boundary, not the fit, set it. "below" is a
# tight binder, "above" a weak one. None when the fit set the value.
Censor = Literal["below", "above"]

# The observed curve's shape. Internal: no Measured/Censored column ships.
CurveClass = Literal["Measured", "Censored High", "Censored Low"]


def gate_index_to_decades(gate_index):
    """Mean gate index → mean log10 fluorescence. Gate label 1 sits at 0 decades.

    Accepts a float, a numpy array or a polars expression.
    """
    return (gate_index - 1.0) * GATE_WIDTH_DECADES


def decades_to_gate_index(decades):
    """The inverse of `gate_index_to_decades`, so fitted and observed curves share one scale."""
    return decades / GATE_WIDTH_DECADES + 1.0


@dataclass(frozen=True)
class GateScale:
    """The fluorescence scale the gates are read on, in decades."""

    floor: float
    top: float

    @property
    def range_decades(self) -> float:
        return self.top - self.floor


def build_gate_scale(max_gate_label: int | None) -> GateScale:
    """The scale implied by the gate labels, at the fixed gate width.

    A single gate label measures no range, so the declared constant stands in rather than dividing
    by zero in the plateau floor.
    """
    floor = 0.0
    top = gate_index_to_decades(float(max_gate_label)) if max_gate_label else floor
    if not top > floor:
        top = floor + GATE_RANGE_DECADES
    return GateScale(floor=floor, top=top)


@dataclass(frozen=True)
class FitResult:
    log10_kd: float | None
    # Standard error of log10 K_D from the fit's own covariance. None when the covariance could not
    # be estimated, which is routine for a fit pinned at a bound.
    log10_kd_se: float | None
    log10_top: float | None  # log10 of the saturated fluorescence B + A
    log10_baseline: float | None
    y_hat: np.ndarray | None  # in decades; convert with decades_to_gate_index before output
    converged: bool
    reason: str | None  # "convergence_failure" on failure; None on success
    r2: float | None = None  # weighted R² on the points the fit saw, in decades


def _failure(reason: str = "convergence_failure") -> FitResult:
    return FitResult(None, None, None, None, None, False, reason)


def _model(c: np.ndarray, log10_kd: float, log10_top: float, log10_base: float) -> np.ndarray:
    base = 10.0**log10_base
    amplitude = 10.0**log10_top - base
    kd = 10.0**log10_kd
    return np.log10(amplitude * c / (c + kd) + base)


def plateau_floor_decades(scale: GateScale, plateau_floor: float) -> float:
    """The log10 fluorescence a saturated clonotype must reach, given the floor as a fraction."""
    return scale.floor + plateau_floor * scale.range_decades


def fit_one(
    c: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    *,
    scale: GateScale,
    plateau_floor: float,
) -> FitResult:
    """Fit one clonotype. Pure numpy; no polars involvement.

    `c` is concentration, `y` the observed mean log10 fluorescence, `w` the reads behind each point.
    The baseline is fitted per clonotype: without absolute cell counts the 0 M mean is not on a
    scale the fit can consume.

    A fit that lands on a bound still converges and is still reported. The caller decides, from
    where K landed relative to the titration, whether the answer is a measurement or a bound.
    """
    c = np.asarray(c, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    positive = c[c > 0]
    if positive.size == 0 or not np.all(np.isfinite(y)):
        return _failure()

    top_lo = plateau_floor_decades(scale, plateau_floor)
    top_hi = max(scale.top, top_lo + _BOUND_GAP)
    base_lo = scale.floor
    base_hi = top_lo - _BOUND_GAP
    if base_hi <= base_lo:
        raise ValueError(
            f"plateau_floor {plateau_floor:.3g} leaves no room between the bottom gate and the "
            f"plateau; raise the floor or widen the gate range."
        )

    log_kd_lo = math.log10(float(positive.min())) - KD_RANGE_DECADES
    log_kd_hi = math.log10(float(positive.max())) + KD_RANGE_DECADES

    y_min, y_max = float(np.min(y)), float(np.max(y))
    top0 = min(max(y_max, top_lo), top_hi)
    base0 = min(max(y_min, base_lo), base_hi)
    half = base0 + (top0 - base0) / 2.0
    above = np.where(y >= half)[0]
    kd0 = float(c[above[0]]) if above.size and c[above[0]] > 0 else float(np.median(positive))
    log_kd0 = min(max(math.log10(kd0), log_kd_lo), log_kd_hi)

    # sigma = 1/sqrt(w); a zero-weight point is treated as very uncertain rather than dropped.
    sigma = np.where(w > 0, 1.0 / np.sqrt(np.maximum(w, 1e-12)), 1e12)

    p0 = [log_kd0, top0, base0]
    lo = [log_kd_lo, top_lo, base_lo]
    hi = [log_kd_hi, top_hi, base_hi]

    try:
        # OptimizeWarning is NOT promoted to an error here. scipy emits it when the covariance
        # cannot be estimated, which is exactly what a fit pinned at a bound does — and a pinned fit
        # is a bound the rating reports as Partial, not a failure.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, pcov = curve_fit(
                _model,
                c,
                y,
                p0=p0,
                sigma=sigma,
                absolute_sigma=False,
                bounds=(lo, hi),
                method="trf",
                maxfev=10_000,
            )
    except (RuntimeError, ValueError):
        return _failure()

    log10_kd = float(popt[0])
    log10_top = float(popt[1])
    log10_base = float(popt[2])
    variance = float(np.asarray(pcov)[0, 0]) if np.ndim(pcov) == 2 else math.nan
    se = math.sqrt(variance) if math.isfinite(variance) and variance >= 0 else None

    y_hat = _model(c, log10_kd, log10_top, log10_base)
    if not np.all(np.isfinite(y_hat)):
        return _failure()

    return FitResult(
        log10_kd=log10_kd,
        log10_kd_se=se,
        log10_top=log10_top,
        log10_baseline=log10_base,
        y_hat=y_hat,
        converged=True,
        reason=None,
        r2=weighted_r2(y, y_hat, w),
    )


def predict(c: np.ndarray, log10_kd: float, log10_top: float, log10_baseline: float) -> np.ndarray:
    """The fitted curve at arbitrary concentrations, in decades."""
    return _model(np.asarray(c, dtype=float), log10_kd, log10_top, log10_baseline)


def rise_test(rise_decades: float | None) -> bool:
    """Did the signal move? Reads the OBSERVED rise, max(F̄) − min(F̄).

    Not the fitted rise: the fitted statistic is pinned at its lower bound for 58 to 76 percent of
    Adams clonotypes, where the bound is set by the gate scale and the baseline rather than by the
    data, so the fitted form flags the tightest binders as non-binders.
    """
    if rise_decades is None or not math.isfinite(rise_decades):
        return False
    return rise_decades >= LF_RISE_THRESHOLD_DECADES


def level_test(plateau_height: float | None, plateau_floor: float) -> bool:
    """For flat curves only. Is the observed top at or above the plateau floor?"""
    if plateau_height is None or not math.isfinite(plateau_height):
        return False
    return plateau_height >= plateau_floor


def classify_curve(
    *,
    rise_decades: float | None,
    plateau_height: float | None,
    plateau_floor: float,
) -> CurveClass:
    """Map one clonotype's observed curve to its class.

    A non-binder is flat and dim. A saturated tight binder is flat and bright. The rise test alone
    calls both non-binders, so the level test — the observed top against the plateau floor — is what
    separates them.
    """
    if rise_test(rise_decades):
        return "Measured"
    if level_test(plateau_height, plateau_floor):
        return "Censored High"
    return "Censored Low"


def pin_kd(
    log10_kd: float,
    *,
    curve_class: CurveClass,
    lo: float,
    hi: float,
) -> tuple[float, Censor | None]:
    """Clamp a K_D to the titrated range and say which wall, if any, set it.

    `lo` and `hi` are log10 of the lowest and highest non-zero concentrations. Returns the reported
    log10 K_D and its censor.
    """
    if curve_class == "Censored High":
        return lo, "below"
    if curve_class == "Censored Low":
        return hi, "above"
    censor: Censor | None = "below" if log10_kd <= lo else "above" if log10_kd >= hi else None
    return min(max(log10_kd, lo), hi), censor
