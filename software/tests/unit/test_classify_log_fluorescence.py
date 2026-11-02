"""The "log_fluorescence" rating of a converged fit, at its boundaries.

Partial  the reported K_D is a bound, or SE(log10 K_D) > KD_SD_GOOD / not finite
Good     K_D inside the titrated range, and SE(log10 K_D) <= KD_SD_GOOD

The pipeline rates gated clonotypes and non-converged fits Failed before this runs; the pipeline
tests cover those.
"""

from __future__ import annotations

import math

import pytest

from classify import classify_log_fluorescence
from constants import KD_SD_GOOD


def _rate(**kwargs):
    base = {"censored": False, "log10_kd_se": 0.1}
    return classify_log_fluorescence(**{**base, **kwargs})


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
