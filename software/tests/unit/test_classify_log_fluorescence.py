"""The "log_fluorescence" rating, at its boundaries.

Failed   the fit did not converge, or too few reads or points to fit
Partial  converged, but the reported K_D is a bound, or SE(log10 K_D) > KD_SD_GOOD / not finite
Good     converged, K_D inside the titrated range, and SE(log10 K_D) <= KD_SD_GOOD
"""

from __future__ import annotations

import math

import pytest

from classify import classify_log_fluorescence
from constants import KD_SD_GOOD


def _rate(**kwargs):
    base = {"converged": True, "censored": False, "log10_kd_se": 0.1}
    return classify_log_fluorescence(**{**base, **kwargs})


class TestFailed:
    def test_a_fit_that_did_not_converge_is_failed(self):
        result = _rate(converged=False)
        assert (result.affinity_class, result.failure_reason) == ("Failed", "convergence_failure")

    @pytest.mark.parametrize("reason", ["insufficient_reads", "insufficient_points"])
    def test_an_upstream_gate_wins_over_everything(self, reason):
        result = _rate(precomputed_reason=reason)
        assert (result.affinity_class, result.failure_reason) == ("Failed", reason)

    def test_a_non_convergence_never_reports_a_fit_reason(self):
        # low_r2 and n_out_of_range cannot fire under this equation.
        assert _rate(converged=False, log10_kd_se=None).failure_reason == "convergence_failure"


class TestGood:
    def test_an_in_range_well_determined_fit_is_good(self):
        result = _rate(log10_kd_se=0.1)
        assert (result.affinity_class, result.failure_reason) == ("Good", None)

    def test_the_standard_error_threshold_is_inclusive(self):
        assert _rate(log10_kd_se=KD_SD_GOOD).affinity_class == "Good"

    def test_a_zero_standard_error_is_good(self):
        assert _rate(log10_kd_se=0.0).affinity_class == "Good"


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
