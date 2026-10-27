"""Behavioral tests for hill_fit.py (R10 Hill kernel, R11 weighted R²)."""

from __future__ import annotations

import math

import numpy as np
import pytest
from helpers import hill_curve
from hypothesis import given
from hypothesis import strategies as st

from constants import DELTA_BIN, DELTA_NO_BIN
from hill_fit import fit_one_clonotype, weighted_r2

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
    # Noiseless round-trip: fit must recover true Kd and n for any curve whose midpoint lies within the
    # concentrations. Guards against reparametrization, sign errors, and wrong Hill algebra.
    @given(log10_kd=st.floats(-1.0, 2.0), n=st.floats(0.5, 2.0), amplitude=st.floats(1.0, 5.0))
    def test_noiseless_roundtrip(self, log10_kd, n, amplitude):
        kd = 10.0**log10_kd
        y = hill_curve(LOG_CONCS, 1.0, math.log(amplitude), kd, n)
        fit = fit_one_clonotype(
            LOG_CONCS, y, np.ones_like(LOG_CONCS), baseline_fixed=1.0, bin_mode=True, max_bin_label=8
        )
        assert fit.converged is True
        assert fit.kd == pytest.approx(kd, rel=1e-4)
        assert fit.n == pytest.approx(n, abs=1e-4)
        assert fit.r2_w == pytest.approx(1.0, abs=1e-6)

    # High-Kd arm: the concentrations barely reach the midpoint, so Kd is only loosely determined.
    def test_noiseless_roundtrip_kd_beyond_the_concentrations(self):
        y = hill_curve(LOG_CONCS, 1.0, math.log(3.0), 1000.0, 1.0)
        fit = fit_one_clonotype(
            LOG_CONCS, y, np.ones_like(LOG_CONCS), baseline_fixed=1.0, bin_mode=True, max_bin_label=8
        )
        assert fit.converged is True
        assert fit.kd == pytest.approx(1000.0, abs=10.0)
        assert fit.n == pytest.approx(1.0, abs=1e-3)

    # No-bin mode: amplitude upper bound is log(0.95); still recoverable.
    def test_no_bin_mode_amplitude_bound(self):
        baseline = 0.01
        amplitude = math.log(0.5)
        x = LOG_CONCS
        y = hill_curve(x, baseline, amplitude, 10.0, 1.0)
        w = np.ones_like(x)

        fit = fit_one_clonotype(x, y, w, baseline_fixed=baseline, bin_mode=False, max_bin_label=None)
        assert fit.converged is True
        assert fit.kd == pytest.approx(10.0, abs=1e-3)

    # 4-parameter fit (B unknown) still recovers baseline within tolerance.
    def test_four_parameter_fit_recovers_baseline(self):
        baseline_true = 1.2
        amplitude = math.log(3.0)
        x = LOG_CONCS
        y = hill_curve(x, baseline_true, amplitude, 10.0, 1.0)
        w = np.ones_like(x)

        fit = fit_one_clonotype(x, y, w, baseline_fixed=None, bin_mode=True, max_bin_label=8)
        assert fit.converged is True
        assert fit.baseline == pytest.approx(baseline_true, abs=1e-3)
        assert fit.kd == pytest.approx(10.0, abs=1e-3)


class TestPinnedHillExponent:
    """`pin_n=True` fixes the Hill exponent at 1, as the Adams 2016 model does.

    `fit_one_clonotype` defaults to `pin_n=False` even though the production default is
    `FitParams.pin_hill_exponent=True`, so a direct call keeps the free-n behaviour.
    """

    @pytest.mark.parametrize(
        "baseline_fixed, bin_mode, max_bin_label, baseline, amplitude",
        [
            pytest.param(1.0, True, 8, 1.0, math.log(3.0), id="bin-baseline-fixed"),
            pytest.param(None, True, 8, 1.0, math.log(3.0), id="bin-baseline-free"),
            pytest.param(0.01, False, None, 0.01, math.log(0.5), id="no-bin-baseline-fixed"),
            pytest.param(None, False, None, 0.01, math.log(0.5), id="no-bin-baseline-free"),
        ],
    )
    def test_pin_n_returns_exactly_one(self, baseline_fixed, bin_mode, max_bin_label, baseline, amplitude):
        # True n = 2.0: a pinned fit must report 1.0, not the data's exponent.
        y = hill_curve(LOG_CONCS, baseline, amplitude, 10.0, 2.0)
        fit = fit_one_clonotype(
            LOG_CONCS,
            y,
            np.full_like(LOG_CONCS, 100.0),
            baseline_fixed,
            bin_mode=bin_mode,
            max_bin_label=max_bin_label,
            pin_n=True,
        )
        assert fit.converged is True
        assert fit.n == 1.0

    def test_default_is_free_n(self):
        y = hill_curve(LOG_CONCS, 1.0, math.log(3.0), 10.0, 2.0)
        fit = fit_one_clonotype(LOG_CONCS, y, np.full_like(LOG_CONCS, 100.0), 1.0, bin_mode=True, max_bin_label=8)
        assert fit.converged is True
        assert fit.n == pytest.approx(2.0, abs=1e-4)


class TestAmplitudePrior:
    """`amp_prior_q` raises the amplitude lower bound to that fraction of the bin range.

    Bin mode with a max_bin_label only. The bound is `1 + q·span` on the saturated mean bin, span
    being `max_bin_label - 1`.
    """

    @pytest.mark.parametrize("baseline_fixed", [1.0, None])
    def test_a_curve_under_the_prior_is_rejected(self, baseline_fixed):
        # True amplitude 1.0 sits well under the prior's floor of 0.65·7 = 4.55.
        y = hill_curve(LOG_CONCS, 1.0, math.log(1.0), 10.0, 1.0)
        w = np.full_like(LOG_CONCS, 100.0)

        free = fit_one_clonotype(LOG_CONCS, y, w, baseline_fixed, bin_mode=True, max_bin_label=8)
        assert free.converged is True
        assert free.top - free.baseline == pytest.approx(1.0, abs=1e-3)

        pinned = fit_one_clonotype(LOG_CONCS, y, w, baseline_fixed, bin_mode=True, max_bin_label=8, amp_prior_q=0.65)
        assert pinned.converged is False
        assert pinned.reason == "insufficient_dynamic_range"

    @pytest.mark.parametrize("baseline_fixed", [1.0, None])
    def test_a_curve_above_the_prior_is_untouched(self, baseline_fixed):
        # True amplitude 6.0 clears the prior's floor of 4.55, so the raised bound never binds.
        y = hill_curve(LOG_CONCS, 1.0, math.log(6.0), 10.0, 1.0)
        w = np.full_like(LOG_CONCS, 100.0)
        fit = fit_one_clonotype(LOG_CONCS, y, w, baseline_fixed, bin_mode=True, max_bin_label=8, amp_prior_q=0.65)
        assert fit.converged is True
        assert fit.top - fit.baseline == pytest.approx(6.0, abs=1e-2)

    # No-bin mode has no bin range, so the floor stays at DELTA_NO_BIN.
    def test_no_bin_mode_ignores_the_prior(self):
        y = hill_curve(LOG_CONCS, 0.01, math.log(0.5), 10.0, 1.0)
        w = np.full_like(LOG_CONCS, 100.0)
        fit = fit_one_clonotype(LOG_CONCS, y, w, 0.01, bin_mode=False, max_bin_label=None, amp_prior_q=0.65)
        assert fit.converged is True
        assert fit.top - fit.baseline == pytest.approx(0.5, abs=1e-6)


class TestHillFitFailureModes:
    # A flat signal reaches the amplitude floor before its weighted R² is ever evaluated: there is
    # no dynamic range for the fit to measure.
    @pytest.mark.parametrize(
        "x, y, baseline_fixed, max_bin_label",
        [
            pytest.param(LOG_CONCS, np.full_like(LOG_CONCS, 1.0), 1.0, 8, id="flat-baseline-fixed"),
            pytest.param(np.array([1e-10, 1e-9, 1e-8, 1e-7, 1e-6]), np.full(5, 2.5), None, 4, id="flat-baseline-free"),
        ],
    )
    def test_flat_signal_has_no_dynamic_range(self, x, y, baseline_fixed, max_bin_label):
        fit = fit_one_clonotype(
            x, y, np.full_like(x, 100.0), baseline_fixed, bin_mode=True, max_bin_label=max_bin_label
        )
        assert not fit.converged
        assert fit.reason == "insufficient_dynamic_range"
        assert fit.kd is None

    # A near-flat ramp reaches the amplitude floor in bin mode too. Admitting these rows instead
    # was measured to cost rank agreement 0.909 -> 0.779 on Adams 2016 replicate 1.
    def test_near_flat_ramp_in_bin_mode_is_rejected(self):
        y = 1.0 + 1e-4 * np.arange(len(LOG_CONCS))
        fit = fit_one_clonotype(LOG_CONCS, y, np.full_like(LOG_CONCS, 100.0), 1.0, bin_mode=True, max_bin_label=8)
        assert fit.converged is False
        assert fit.reason == "insufficient_dynamic_range"
        assert fit.kd is None

    # The same ramp in no-bin mode, where the floor is an absolute frequency.
    def test_near_flat_ramp_in_no_bin_mode_is_rejected(self):
        y = 1e-4 * np.arange(len(LOG_CONCS))
        fit = fit_one_clonotype(LOG_CONCS, y, np.full_like(LOG_CONCS, 100.0), 0.0, bin_mode=False, max_bin_label=None)
        assert fit.converged is False
        assert fit.reason == "insufficient_dynamic_range"
        assert fit.kd is None


class TestDeltaDynamicRangeGate:
    """R10: `top − baseline` must be ≥ mode-specific δ (δ_bin = 0.5, δ_no_bin = 0.05).

    Bin mode applies δ as a bound on the optimiser only: a curve below δ comes back converged with
    its amplitude pinned at δ. No-bin mode keeps the post-fit rejection, with reason
    `insufficient_dynamic_range`.
    """

    @pytest.mark.parametrize("amplitude", [0.05, 0.1, 0.3])
    def test_bin_mode_rejects_a_sub_delta_curve(self, amplitude):
        y = hill_curve(LOG_CONCS, 1.0, math.log(amplitude), kd=10.0, n=1.0)
        fit = fit_one_clonotype(
            LOG_CONCS, y, np.full_like(LOG_CONCS, 100.0), baseline_fixed=1.0, bin_mode=True, max_bin_label=8
        )
        assert fit.converged is False
        assert fit.reason == "insufficient_dynamic_range"
        assert fit.kd is None

    @pytest.mark.parametrize("amplitude", [0.6, 2.0])
    def test_bin_mode_recovers_a_curve_above_delta(self, amplitude):
        y = hill_curve(LOG_CONCS, 1.0, math.log(amplitude), kd=10.0, n=1.0)
        fit = fit_one_clonotype(
            LOG_CONCS, y, np.full_like(LOG_CONCS, 100.0), baseline_fixed=1.0, bin_mode=True, max_bin_label=8
        )
        assert fit.converged is True
        assert fit.top - fit.baseline == pytest.approx(amplitude, abs=1e-4)

    @pytest.mark.parametrize("amplitude", [0.01, 0.03])
    def test_no_bin_mode_rejects_a_sub_delta_curve(self, amplitude):
        y = hill_curve(LOG_CONCS, 0.01, math.log(amplitude), kd=10.0, n=1.0)
        fit = fit_one_clonotype(
            LOG_CONCS, y, np.full_like(LOG_CONCS, 100.0), baseline_fixed=0.01, bin_mode=False, max_bin_label=None
        )
        assert fit.converged is False
        assert fit.reason == "insufficient_dynamic_range"
        assert fit.kd is None

    @pytest.mark.parametrize("amplitude", [0.08, 0.5])
    def test_no_bin_mode_recovers_a_curve_above_delta(self, amplitude):
        y = hill_curve(LOG_CONCS, 0.01, math.log(amplitude), kd=10.0, n=1.0)
        fit = fit_one_clonotype(
            LOG_CONCS, y, np.full_like(LOG_CONCS, 100.0), baseline_fixed=0.01, bin_mode=False, max_bin_label=None
        )
        assert fit.converged is True
        assert fit.top - fit.baseline == pytest.approx(amplitude, abs=1e-4)


class TestOptimizeWarningBecomesFailure:
    """Contract: a scipy OptimizeWarning during curve_fit means the fit is unreliable —
    treat it as convergence_failure rather than returning a point estimate whose
    covariance couldn't be estimated.

    We mock the external boundary (scipy.optimize.curve_fit inside hill_fit) so the
    test is deterministic across scipy versions and test machines, where natural
    triggers (singular Jacobians) vary.
    """

    def test_optimize_warning_returns_failure(self, monkeypatch):
        import warnings as _warnings

        from scipy.optimize import OptimizeWarning

        import hill_fit

        def fake_curve_fit(*_args: object, **_kwargs: object) -> tuple[np.ndarray, None]:
            _warnings.warn(
                "Covariance of the parameters could not be estimated",
                OptimizeWarning,
                stacklevel=2,
            )
            # Realistic popt shape for baseline-fixed mode: [log_kd, n, amplitude].
            return np.array([math.log(10.0), 1.0, math.log(2.0)]), None

        monkeypatch.setattr(hill_fit, "curve_fit", fake_curve_fit)

        x = np.array([1e-9, 1e-8, 1e-7, 1e-6])
        y = np.array([1.0, 1.5, 2.5, 3.0])
        w = np.ones_like(x) * 100.0
        fit = fit_one_clonotype(x, y, w, baseline_fixed=1.0, bin_mode=True, max_bin_label=8)

        assert fit.converged is False
        assert fit.reason == "convergence_failure"


class TestAmplitudeBoundsRespected:
    """R10 reparametrization: every CONVERGED fit must have top − baseline ≥ δ_mode.

    Belt-and-braces invariant — guards against future refactors of the amp_lo
    solver bound drifting back below the spec minimum. Complements TestDeltaDynamicRangeGate.
    """

    @pytest.mark.parametrize("bin_mode, delta_expected", [(True, DELTA_BIN), (False, DELTA_NO_BIN)])
    @pytest.mark.parametrize("true_amplitude_fraction", [0.1, 1.0])
    def test_converged_top_minus_baseline_ge_delta(self, bin_mode, delta_expected, true_amplitude_fraction):
        baseline = 1.0 if bin_mode else 0.01
        max_bin_label = 8 if bin_mode else None
        amplitude = math.log((1.5 if bin_mode else 0.3) * true_amplitude_fraction)
        x = LOG_CONCS
        y = hill_curve(x, baseline, amplitude, kd=10.0, n=1.0)
        w = np.ones_like(x) * 100.0

        fit = fit_one_clonotype(x, y, w, baseline_fixed=baseline, bin_mode=bin_mode, max_bin_label=max_bin_label)
        if not fit.converged:
            assert fit.reason == "insufficient_dynamic_range"
            return
        assert fit.top is not None and fit.baseline is not None
        assert (fit.top - fit.baseline) >= delta_expected - 1e-9, (
            f"top − baseline = {fit.top - fit.baseline:.4f} violates δ = {delta_expected}"
        )
