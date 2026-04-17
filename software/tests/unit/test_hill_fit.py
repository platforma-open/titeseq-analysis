"""Behavioral tests for hill_fit.py (R10 Hill kernel, R11 weighted R²)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from hill_fit import fit_one_clonotype, weighted_r2


def hill_truth(x, baseline, amplitude, kd, n):
    top = baseline + math.exp(amplitude)
    return baseline + (top - baseline) * (x**n) / (kd**n + x**n)


LOG_CONCS = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0])


class TestWeightedR2:
    # Pinned values for weighted R². Guards against clamping / ordering bugs.
    @pytest.mark.parametrize(
        "y, y_hat, w, expected",
        [
            ([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 1.0, 1.0], 1.0),  # perfect fit
            ([1.0, 2.0, 3.0], [2.0, 2.0, 2.0], [1.0, 1.0, 1.0], 0.0),  # fit = constant mean
            ([1.0, 2.0, 3.0], [3.0, 2.0, 1.0], [1.0, 1.0, 1.0], -3.0),  # anti-fit, NOT clamped
            # Zero-weight point: excluded from both num and weighted mean.
            # ŷ matches y on weighted points → R² = 1.0 regardless of y_hat[2].
            ([1.0, 2.0, 3.0], [1.0, 2.0, 100.0], [1.0, 1.0, 0.0], 1.0),
            # Partial-weight:
            # ȳ_w = (1·1 + 1·2 + 0.5·3) / 2.5 = 1.8
            # SSR = 0; SST = 1·0.64 + 1·0.04 + 0.5·1.44 = 1.40 ≠ 0 → R² = 1.0
            ([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 1.0, 0.5], 1.0),
        ],
    )
    def test_weighted_r2_pinned(self, y, y_hat, w, expected):
        assert weighted_r2(np.asarray(y), np.asarray(y_hat), np.asarray(w)) == pytest.approx(expected)

    # Degenerate inputs (Σw = 0, zero-variance y under weighting) return nan sentinel.
    @pytest.mark.parametrize(
        "y, y_hat, w",
        [
            ([1.0, 2.0], [1.0, 2.0], [0.0, 0.0]),  # Σw = 0
            ([2.0, 2.0, 2.0], [2.0, 2.0, 2.0], [1.0, 1.0, 1.0]),  # zero-variance y
        ],
    )
    def test_degenerate_returns_nan(self, y, y_hat, w):
        assert math.isnan(weighted_r2(np.array(y), np.array(y_hat), np.array(w)))


class TestHillFitRoundtrip:
    # Noiseless round-trip: fit must recover true K_D within tight tolerance.
    # Guards against reparametrization, sign errors, and wrong Hill algebra.
    @pytest.mark.parametrize(
        "true_kd, true_n, abs_err_kd, abs_err_n",
        [
            (10.0, 1.0, 1e-4, 1e-4),
            (10.0, 2.0, 1e-4, 1e-4),
            (10.0, 0.5, 1e-3, 1e-4),
            (0.1, 1.0, 1e-4, 1e-4),
            (1000.0, 1.0, 10.0, 1e-3),  # high-KD arm: concentrations barely reach mid-point
        ],
    )
    def test_noiseless_roundtrip(self, true_kd, true_n, abs_err_kd, abs_err_n):
        baseline = 1.0
        amplitude = math.log(3.0)
        x = LOG_CONCS
        y = hill_truth(x, baseline, amplitude, true_kd, true_n)
        w = np.ones_like(x)

        fit = fit_one_clonotype(x, y, w, baseline_fixed=baseline, bin_mode=True, max_bin_label=8)
        assert fit.converged is True
        assert fit.kd == pytest.approx(true_kd, abs=abs_err_kd, rel=1e-3)
        assert fit.n == pytest.approx(true_n, abs=abs_err_n)
        assert fit.r2_w == pytest.approx(1.0, abs=1e-6)

    # No-bin mode: amplitude upper bound is log(0.95); still recoverable.
    def test_no_bin_mode_amplitude_bound(self):
        baseline = 0.01
        amplitude = math.log(0.5)
        x = LOG_CONCS
        y = hill_truth(x, baseline, amplitude, 10.0, 1.0)
        w = np.ones_like(x)

        fit = fit_one_clonotype(x, y, w, baseline_fixed=baseline, bin_mode=False, max_bin_label=None)
        assert fit.converged is True
        assert fit.kd == pytest.approx(10.0, abs=1e-3)

    # 4-parameter fit (B unknown) still recovers baseline within tolerance.
    def test_four_parameter_fit_recovers_baseline(self):
        baseline_true = 1.2
        amplitude = math.log(3.0)
        x = LOG_CONCS
        y = hill_truth(x, baseline_true, amplitude, 10.0, 1.0)
        w = np.ones_like(x)

        fit = fit_one_clonotype(x, y, w, baseline_fixed=None, bin_mode=True, max_bin_label=8)
        assert fit.converged is True
        assert fit.baseline == pytest.approx(baseline_true, abs=1e-3)
        assert fit.kd == pytest.approx(10.0, abs=1e-3)


class TestHillFitFailureModes:
    # Flat curve: top − baseline < δ → convergence_failure.
    # This specifically exercises the `top - baseline < δ` branch of R10.
    def test_flat_curve_yields_convergence_failure(self):
        x = LOG_CONCS
        y = np.full_like(x, 1.0)  # absolutely flat signal
        w = np.ones_like(x)
        fit = fit_one_clonotype(x, y, w, baseline_fixed=1.0, bin_mode=True, max_bin_label=8)
        assert fit.converged is False
        assert fit.reason == "convergence_failure"

    # Tiny dynamic range (< δ) with baseline fixed → reject with convergence_failure.
    def test_tiny_amplitude_below_delta_fails(self):
        x = LOG_CONCS
        # Amplitude well below δ=0.05 in absolute terms
        y = np.full_like(x, 1.0) + 1e-4 * np.arange(len(x))
        w = np.ones_like(x)
        fit = fit_one_clonotype(x, y, w, baseline_fixed=1.0, bin_mode=True, max_bin_label=8)
        assert fit.converged is False
        assert fit.reason == "convergence_failure"


class TestDeltaDynamicRangeGate:
    """R10: `top − baseline` must be ≥ mode-specific δ.

    Spec: δ_bin = 0.5, δ_no_bin = 0.05. Fits whose true dynamic range falls
    below the mode's δ must be rejected as convergence_failure; fits above
    must converge.
    """

    @pytest.mark.parametrize(
        "bin_mode, amplitude, should_converge",
        [
            # Bin mode, δ = 0.5:
            (True, math.log(0.3), False),  # 0.3 < 0.5 → reject (current code ACCEPTS; catches bug)
            (True, math.log(0.6), True),  # 0.6 ≥ 0.5 → accept
            (True, math.log(2.0), True),  # healthy signal
            # No-bin mode, δ = 0.05:
            (False, math.log(0.03), False),  # 0.03 < 0.05 → reject
            (False, math.log(0.08), True),  # 0.08 ≥ 0.05 → accept
        ],
    )
    def test_mode_specific_delta_gate(self, bin_mode, amplitude, should_converge):
        baseline = 1.0 if bin_mode else 0.01
        max_bin_label = 8 if bin_mode else None
        x = LOG_CONCS
        y = hill_truth(x, baseline, amplitude, kd=10.0, n=1.0)
        w = np.ones_like(x) * 100.0

        fit = fit_one_clonotype(
            x, y, w, baseline_fixed=baseline, bin_mode=bin_mode, max_bin_label=max_bin_label
        )
        assert fit.converged is should_converge, (
            f"{'bin' if bin_mode else 'no-bin'} amp={math.exp(amplitude):.3f}: "
            f"expected converged={should_converge}, got {fit.converged}"
        )
        if not should_converge:
            assert fit.reason == "convergence_failure"


class TestAmplitudeBoundsRespected:
    """R10 reparametrization: every CONVERGED fit must have top − baseline ≥ δ_mode.

    Belt-and-braces invariant — guards against future refactors of the amp_lo
    solver bound drifting back below the spec minimum. Complements TestDeltaDynamicRangeGate.
    """

    @pytest.mark.parametrize("bin_mode, delta_expected", [(True, 0.5), (False, 0.05)])
    def test_converged_top_minus_baseline_ge_delta(self, bin_mode, delta_expected):
        baseline = 1.0 if bin_mode else 0.01
        max_bin_label = 8 if bin_mode else None
        amplitude = math.log(1.5) if bin_mode else math.log(0.3)
        x = LOG_CONCS
        y = hill_truth(x, baseline, amplitude, kd=10.0, n=1.0)
        w = np.ones_like(x) * 100.0

        fit = fit_one_clonotype(
            x, y, w, baseline_fixed=baseline, bin_mode=bin_mode, max_bin_label=max_bin_label
        )
        assert fit.converged is True
        assert (fit.top - fit.baseline) >= delta_expected - 1e-9, (
            f"top − baseline = {fit.top - fit.baseline:.4f} violates δ = {delta_expected}"
        )
