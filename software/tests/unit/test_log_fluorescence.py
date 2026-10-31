"""The "log_fluorescence" equation: the gate scale, parameter recovery, the plateau floor, the pin."""

from __future__ import annotations

import math

import numpy as np
import pytest

from constants import GATE_RANGE_DECADES, GATE_WIDTH_DECADES, LF_RISE_THRESHOLD_DECADES
from log_fluorescence import (
    GateScale,
    build_gate_scale,
    classify_curve,
    decades_to_gate_index,
    fit_one,
    gate_index_to_decades,
    pin_kd,
    plateau_floor_decades,
    predict,
)

CONCS = np.array([1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6])
# Four gates, 0.65 decades apart: the scale the block assumes.
SCALE = GateScale(floor=0.0, top=1.95)
FLOOR = 2.0 / 3.0
LOG_LO = math.log10(CONCS.min())
LOG_HI = math.log10(CONCS.max())


def _curve(kd: float, *, baseline: float = 2.0, top: float = 10 ** (FLOOR * 1.95 + 0.3)) -> np.ndarray:
    """log10(A*c/(c+K) + B) — the model the block fits, written out so the test does not reuse src."""
    return np.log10((top - baseline) * CONCS / (CONCS + kd) + baseline)


def _fit(y: np.ndarray, *, w: np.ndarray | None = None, plateau_floor: float = FLOOR):
    weights = np.full(CONCS.size, 1000.0) if w is None else w
    return fit_one(CONCS, y, weights, scale=SCALE, plateau_floor=plateau_floor)


class TestGateScale:
    def test_gate_index_and_decades_are_inverses(self):
        assert gate_index_to_decades(1.0) == pytest.approx(0.0)
        assert gate_index_to_decades(4.0) == pytest.approx(3 * GATE_WIDTH_DECADES)
        assert decades_to_gate_index(gate_index_to_decades(2.7)) == pytest.approx(2.7)

    def test_four_gates_span_three_gate_widths(self):
        scale = build_gate_scale(4)
        assert scale.floor == pytest.approx(0.0)
        assert scale.top == pytest.approx(1.95)
        assert scale.range_decades == pytest.approx(1.95)

    @pytest.mark.parametrize("max_gate", [None, 0, 1])
    def test_a_scale_that_measures_nothing_falls_back_to_the_declared_range(self, max_gate):
        assert build_gate_scale(max_gate).range_decades == pytest.approx(GATE_RANGE_DECADES)

    def test_plateau_floor_decades_is_a_fraction_of_the_gate_range(self):
        assert plateau_floor_decades(SCALE, 2.0 / 3.0) == pytest.approx(1.3)
        assert plateau_floor_decades(GateScale(1.0, 3.0), 0.5) == pytest.approx(2.0)


class TestFit:
    @pytest.mark.parametrize("log10_kd", [-9.0, -8.0, -7.0, -6.5])
    def test_recovers_a_known_kd_from_a_noiseless_curve(self, log10_kd):
        fit = _fit(_curve(10.0**log10_kd))
        assert fit.converged
        assert fit.log10_kd == pytest.approx(log10_kd, abs=0.05)

    def test_recovers_the_baseline_it_was_never_given(self):
        fit = _fit(_curve(1e-8, baseline=2.5))
        assert fit.converged and fit.log10_baseline is not None
        assert fit.log10_kd == pytest.approx(-8.0, abs=0.05)
        assert fit.log10_baseline == pytest.approx(math.log10(2.5), abs=0.05)

    def test_the_exponent_is_fixed_at_one(self):
        # A curve with a steeper exponent cannot be matched, so the fit's own residuals carry it.
        y = np.log10(8.0 * CONCS**2 / (CONCS**2 + (1e-8) ** 2) + 2.0)
        fit = _fit(y)
        assert fit.converged
        assert fit.r2 is not None and fit.r2 < 1.0

    def test_the_fitted_curve_reproduces_the_points(self):
        y = _curve(1e-8)
        fit = _fit(y)
        assert fit.y_hat is not None
        assert fit.log10_kd is not None and fit.log10_top is not None and fit.log10_baseline is not None
        assert np.allclose(fit.y_hat, y, atol=0.02)
        assert np.allclose(predict(CONCS, fit.log10_kd, fit.log10_top, fit.log10_baseline), fit.y_hat)
        assert fit.r2 == pytest.approx(1.0, abs=1e-3)

    def test_a_saturated_curve_is_reported_below_the_titration_not_as_a_failure(self):
        # Saturated everywhere: K_D is below the titration. The fit still converges — the bound is
        # what makes the row Partial, not a refusal.
        fit = _fit(_curve(1e-14))
        assert fit.converged and fit.log10_kd is not None
        assert fit.log10_kd < LOG_LO

    def test_a_flat_dim_curve_is_pushed_above_the_titration_by_the_plateau_floor(self):
        # A flat non-binder has no rise to fit. Without the floor the optimiser would place a tiny
        # amplitude anywhere; the floor forces the only consistent answer, K_D above the range.
        fit = _fit(np.full(CONCS.size, 0.05))
        assert fit.converged and fit.log10_kd is not None
        assert fit.log10_kd > LOG_HI

    def test_a_lower_plateau_floor_admits_a_lower_ceiling_as_a_real_binder(self):
        # A clonotype that saturates at 45 percent of the gate range: a non-binder at the screening
        # floor, a measured binder at the lead-optimisation floor.
        y = _curve(1e-8, top=10 ** (0.45 * 1.95))
        lenient = _fit(y, plateau_floor=0.4)
        strict = _fit(y, plateau_floor=FLOOR)
        assert lenient.log10_kd is not None and strict.log10_kd is not None
        assert lenient.log10_kd == pytest.approx(-8.0, abs=0.1)
        assert strict.log10_kd > lenient.log10_kd

    def test_zero_weight_points_do_not_break_the_fit(self):
        w = np.full(CONCS.size, 1000.0)
        w[0] = 0.0
        fit = _fit(_curve(1e-8), w=w)
        assert fit.converged and fit.log10_kd is not None
        assert fit.log10_kd == pytest.approx(-8.0, abs=0.05)

    def test_a_standard_error_comes_back_for_a_well_determined_fit(self):
        noise = np.array([0.01, -0.01, 0.02, -0.02, 0.0, 0.01, -0.01, 0.0, 0.01])
        fit = _fit(_curve(1e-8) + noise)
        assert fit.log10_kd_se is not None
        assert 0.0 < fit.log10_kd_se < 0.5

    def test_a_plateau_floor_leaving_no_baseline_room_is_refused_loudly(self):
        with pytest.raises(ValueError, match="plateau_floor"):
            _fit(_curve(1e-8), plateau_floor=1e-4)

    def test_no_positive_concentration_is_a_failure_not_a_crash(self):
        fit = fit_one(np.zeros(4), np.zeros(4), np.ones(4), scale=SCALE, plateau_floor=FLOOR)
        assert not fit.converged
        assert fit.reason == "convergence_failure"
        assert fit.log10_kd is None and fit.r2 is None

    def test_a_non_finite_signal_is_a_failure_not_a_crash(self):
        y = _curve(1e-8)
        y[2] = math.nan
        fit = _fit(y)
        assert not fit.converged
        assert fit.reason == "convergence_failure"


class TestCurveClass:
    @pytest.mark.parametrize("rise", [LF_RISE_THRESHOLD_DECADES, LF_RISE_THRESHOLD_DECADES + 0.01, 2.0])
    def test_a_curve_that_rose_far_enough_is_measured(self, rise):
        assert classify_curve(rise_decades=rise, plateau_height=0.9, plateau_floor=FLOOR) == "Measured"

    def test_the_rise_threshold_is_inclusive_from_below(self):
        assert (
            classify_curve(
                rise_decades=LF_RISE_THRESHOLD_DECADES - 1e-9,
                plateau_height=0.9,
                plateau_floor=FLOOR,
            )
            == "Censored High"
        )

    def test_a_flat_bright_curve_is_censored_high(self):
        assert classify_curve(rise_decades=0.0, plateau_height=FLOOR, plateau_floor=FLOOR) == "Censored High"

    def test_a_flat_dim_curve_is_censored_low(self):
        assert classify_curve(rise_decades=0.0, plateau_height=FLOOR - 1e-9, plateau_floor=FLOOR) == "Censored Low"

    @pytest.mark.parametrize("bad", [None, math.nan, math.inf])
    def test_an_unusable_statistic_never_claims_a_measurement(self, bad):
        assert classify_curve(rise_decades=bad, plateau_height=bad, plateau_floor=FLOOR) == "Censored Low"


class TestPin:
    def test_a_measurement_inside_the_range_is_untouched(self):
        value, censor = pin_kd(-8.0, curve_class="Measured", lo=LOG_LO, hi=LOG_HI)
        assert value == pytest.approx(-8.0)
        assert censor is None

    def test_a_measurement_below_the_range_is_pinned_to_the_low_wall(self):
        value, censor = pin_kd(-14.0, curve_class="Measured", lo=LOG_LO, hi=LOG_HI)
        assert value == pytest.approx(LOG_LO)
        assert censor == "below"

    def test_a_measurement_above_the_range_is_pinned_to_the_high_wall(self):
        value, censor = pin_kd(-2.0, curve_class="Measured", lo=LOG_LO, hi=LOG_HI)
        assert value == pytest.approx(LOG_HI)
        assert censor == "above"

    @pytest.mark.parametrize("edge,censor", [(LOG_LO, "below"), (LOG_HI, "above")])
    def test_a_value_exactly_on_a_wall_counts_as_set_by_that_wall(self, edge, censor):
        value, got = pin_kd(edge, curve_class="Measured", lo=LOG_LO, hi=LOG_HI)
        assert value == pytest.approx(edge)
        assert got == censor

    def test_censored_high_reports_the_lowest_concentration(self):
        assert pin_kd(0.0, curve_class="Censored High", lo=LOG_LO, hi=LOG_HI) == (LOG_LO, "below")

    def test_censored_low_reports_the_highest_concentration(self):
        assert pin_kd(0.0, curve_class="Censored Low", lo=LOG_LO, hi=LOG_HI) == (LOG_HI, "above")
