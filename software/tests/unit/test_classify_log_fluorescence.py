"""The "log_fluorescence" rating of a converged fit, at its boundaries.

Non-binder  a flat, dim curve, its K_D pinned at the weak wall
Partial     the reported K_D is a bound
Failed      low_r2: the R² gate is on and an in-range fit's R² is below it or undefined
Partial     SE(log10 K_D) > KD_SD_GOOD / not finite
Good        K_D inside the titrated range, and SE(log10 K_D) <= KD_SD_GOOD

The pipeline rates gated clonotypes and non-converged fits Failed before this runs; the pipeline
tests cover those.
"""

from __future__ import annotations

import math

import pytest

from classify import classify_log_fluorescence
from constants import KD_SD_GOOD


def _rate(**kwargs):
    base = {"non_binder": False, "censored": False, "log10_kd_se": 0.1, "r2": 0.95, "r2_gate_min": None}
    return classify_log_fluorescence(**{**base, **kwargs})


class TestGood:
    def test_an_in_range_well_determined_fit_is_good(self):
        result = _rate(log10_kd_se=0.1)
        assert (result.affinity_class, result.failure_reason) == ("Good", None)

    def test_the_standard_error_threshold_is_inclusive(self):
        assert _rate(log10_kd_se=KD_SD_GOOD).affinity_class == "Good"

    def test_a_zero_standard_error_is_good(self):
        assert _rate(log10_kd_se=0.0).affinity_class == "Good"


class TestNonBinder:
    def test_a_flat_dim_curve_is_a_non_binder(self):
        result = _rate(non_binder=True, censored=True)
        assert (result.affinity_class, result.failure_reason) == ("Non-binder", None)

    def test_the_r2_gate_does_not_reach_a_non_binder(self):
        assert _rate(non_binder=True, censored=True, r2=-1.0, r2_gate_min=0.8).affinity_class == "Non-binder"


class TestPartial:
    def test_a_bound_is_partial_however_tight_the_fit_looked(self):
        result = _rate(censored=True, log10_kd_se=0.0)
        assert (result.affinity_class, result.failure_reason) == ("Partial", None)

    def test_just_above_the_standard_error_threshold_is_partial(self):
        assert _rate(log10_kd_se=KD_SD_GOOD + 1e-9).affinity_class == "Partial"

    @pytest.mark.parametrize("se", [None, math.nan, math.inf])
    def test_a_standard_error_that_is_not_a_finite_number_is_partial(self, se):
        result = _rate(log10_kd_se=se)
        assert (result.affinity_class, result.failure_reason) == ("Partial", None)

    def test_partial_never_carries_a_failure_reason(self):
        assert _rate(censored=True).failure_reason is None
        assert _rate(log10_kd_se=None).failure_reason is None


class TestR2Gate:
    def test_off_by_default_a_low_r2_is_still_good(self):
        assert _rate(r2=-0.5).affinity_class == "Good"

    def test_an_in_range_fit_below_the_gate_fails_low_r2(self):
        result = _rate(r2=0.79, r2_gate_min=0.8)
        assert (result.affinity_class, result.failure_reason) == ("Failed", "low_r2")

    def test_the_gate_is_inclusive(self):
        assert _rate(r2=0.8, r2_gate_min=0.8).affinity_class == "Good"

    @pytest.mark.parametrize("r2", [None, math.nan])
    def test_an_undefined_r2_fails_when_the_gate_is_on(self, r2):
        assert _rate(r2=r2, r2_gate_min=0.8).failure_reason == "low_r2"

    def test_a_bound_is_not_gated(self):
        # Phillips 2021 filters in-range fits only.
        assert _rate(censored=True, r2=0.1, r2_gate_min=0.8).affinity_class == "Partial"

    def test_the_gate_runs_before_the_standard_error(self):
        assert _rate(r2=0.5, r2_gate_min=0.8, log10_kd_se=5.0).affinity_class == "Failed"


class TestSeGate:
    def test_off_by_default_a_wide_error_is_partial(self):
        assert _rate(log10_kd_se=3.0).affinity_class == "Partial"

    def test_an_in_range_fit_above_the_limit_fails_high_kd_se(self):
        result = _rate(log10_kd_se=1.01, kd_se_gate_max=1.0)
        assert (result.affinity_class, result.failure_reason) == ("Failed", "high_kd_se")

    def test_the_limit_is_inclusive(self):
        assert _rate(log10_kd_se=1.0, kd_se_gate_max=1.0).affinity_class == "Partial"

    @pytest.mark.parametrize("se", [None, math.nan])
    def test_an_undefined_se_fails_when_the_limit_is_set(self, se):
        assert _rate(log10_kd_se=se, kd_se_gate_max=1.0).failure_reason == "high_kd_se"

    def test_a_bound_is_not_gated(self):
        assert _rate(censored=True, log10_kd_se=None, kd_se_gate_max=1.0).affinity_class == "Partial"

    def test_the_r2_gate_reads_first(self):
        result = _rate(r2=0.1, r2_gate_min=0.8, log10_kd_se=5.0, kd_se_gate_max=1.0)
        assert result.failure_reason == "low_r2"
