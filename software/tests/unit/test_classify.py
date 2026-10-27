"""Behavioral tests for classify.py (R12 truth table + boundaries)."""

from __future__ import annotations

import math
from typing import get_args

import pytest

from classify import classify
from constants import DEFAULT_PARAMS, AffinityClass, FailureReason, FitParams


# Truth table: a converged fit is never Failed. R² and n set Good vs Partial; neither discards it.
@pytest.mark.parametrize(
    "r2, n, converged, expected_class, expected_reason",
    [
        # R² >= Good, n in range → Good
        (0.90, 1.0, True, "Good", None),
        (0.80, 1.0, True, "Good", None),  # boundary: R² == r2_threshold_good (inclusive)
        # R² >= Good, n out of range → Partial
        (0.95, 0.4, True, "Partial", None),
        (0.95, 2.1, True, "Partial", None),
        (0.95, 0.5, True, "Good", None),  # boundary: n == n_min (inclusive)
        (0.95, 2.0, True, "Good", None),  # boundary: n == n_max (inclusive)
        # R² below Good, n in range → Partial
        (0.79, 1.0, True, "Partial", None),  # boundary: just under r2_threshold_good
        (0.50, 1.0, True, "Partial", None),
        # R² below Good, n out of range → Partial (no n_out_of_range outcome exists)
        (0.70, 0.3, True, "Partial", None),
        (0.70, 2.5, True, "Partial", None),
        # Terrible R² on a converged fit is still reported (no low_r2 outcome exists)
        (0.49, 1.0, True, "Partial", None),
        (0.0, 1.0, True, "Partial", None),
        (-0.5, 1.0, True, "Partial", None),
        (-1e6, 1.0, True, "Partial", None),
        (-0.5, 5.0, True, "Partial", None),  # terrible R² and n far out of range
        # A non-finite R² (flat signal, no variance to explain) clears no threshold → Partial
        (float("nan"), 1.0, True, "Partial", None),
        (float("inf"), 1.0, True, "Partial", None),
        (float("-inf"), 1.0, True, "Partial", None),
        # Did not converge → Failed (convergence_failure), regardless of R²/n
        (0.0, 0.0, False, "Failed", "convergence_failure"),
        (0.99, 1.0, False, "Failed", "convergence_failure"),
        (None, None, False, "Failed", "convergence_failure"),
        # Converged, but the fit left R² or n undefined → Failed (convergence_failure)
        (None, 1.0, True, "Failed", "convergence_failure"),
        (0.95, None, True, "Failed", "convergence_failure"),
    ],
)
def test_classification_truth_table(r2, n, converged, expected_class, expected_reason):
    result = classify(r2, n, converged, DEFAULT_PARAMS)
    assert result.affinity_class == expected_class
    assert result.failure_reason == expected_reason


# A tight binder whose curve the fitter explains badly must still reach the table with its K_D.
# The adequacy gate that discarded these cost 215/209/217 of the 250 tightest binders on the three
# Adams 2016 replicates.
@pytest.mark.parametrize("r2", [0.49, 0.0, -0.5, -12.0, float("nan")])
def test_converged_fit_is_never_failed(r2):
    result = classify(r2, 1.0, True, DEFAULT_PARAMS)
    assert result.affinity_class != "Failed"
    assert result.failure_reason is None


# A flat signal reaches classify with an undefined or near-zero R²: never a confident binder.
@pytest.mark.parametrize("r2", [float("nan"), 0.0, -0.5])
def test_flat_signal_is_not_reported_as_a_confident_binder(r2):
    assert classify(r2, 1.0, True, DEFAULT_PARAMS).affinity_class != "Good"


# Non-default thresholds (e.g. multimeric-antigen scenario) exercise the same truth table.
# Guards against hardcoded constants leaking into classify().
def test_classification_with_non_default_thresholds():
    strict = FitParams(r2_threshold_good=0.95, r2_threshold_failed=0.7, n_min=0.8, n_max=3.0)
    # R² = 0.93, n=1.5: in Partial zone under strict thresholds (would be Good with defaults).
    result = classify(0.93, 1.5, True, strict)
    assert result.affinity_class == "Partial"
    assert result.failure_reason is None

    # R² = 0.96, n=0.5: n_min is 0.8 so n is out of range → Partial.
    result = classify(0.96, 0.5, True, strict)
    assert result.affinity_class == "Partial"

    # n = 2.5 is out of range under the defaults and in range under `strict`.
    assert classify(0.96, 2.5, True, strict).affinity_class == "Good"
    assert classify(0.96, 2.5, True, DEFAULT_PARAMS).affinity_class == "Partial"


# Precomputed upstream gates (insufficient_*, non_monotonic_signal) win over fit-based reasons.
@pytest.mark.parametrize(
    "precomputed",
    [
        "insufficient_reads",
        "insufficient_points",
        "non_monotonic_signal",
        "insufficient_dynamic_range",
        "convergence_failure",
    ],
)
def test_precomputed_reason_wins(precomputed):
    result = classify(0.99, 1.0, True, DEFAULT_PARAMS, precomputed_reason=precomputed)
    assert result.affinity_class == "Failed"
    assert result.failure_reason == precomputed


# A refusal upstream of the fit keeps its own reason rather than being relabelled.
def test_precomputed_reason_survives_a_non_converged_fit():
    result = classify(None, None, False, DEFAULT_PARAMS, precomputed_reason="insufficient_dynamic_range")
    assert result.affinity_class == "Failed"
    assert result.failure_reason == "insufficient_dynamic_range"


class TestConstants:
    # Every reason is a refusal to fit at all; "low_r2", "n_out_of_range" and "kd_above_range" are gone.
    def test_failure_reason_vocabulary(self):
        assert get_args(FailureReason) == (
            "insufficient_reads",
            "insufficient_points",
            "non_monotonic_signal",
            "convergence_failure",
            "insufficient_dynamic_range",
        )

    def test_affinity_class_vocabulary(self):
        assert get_args(AffinityClass) == ("Good", "Partial", "CensoredBelow", "Failed")

    def test_fit_params_defaults(self):
        assert DEFAULT_PARAMS.hook_detector_enabled is False
        assert DEFAULT_PARAMS.pin_hill_exponent is True
        assert DEFAULT_PARAMS.amplitude_prior_q is None

    @pytest.mark.parametrize("removed", ["allow_computed_baseline", "allow_uncorrected"])
    def test_removed_fit_params_fields_are_gone(self, removed):
        assert not hasattr(DEFAULT_PARAMS, removed)


# classify() emits no fit-based reason for any converged input, finite or not.
@pytest.mark.parametrize("r2", [1.0, 0.8, 0.5, 0.0, -3.0, float("nan"), float("inf")])
@pytest.mark.parametrize("n", [0.1, 0.5, 1.0, 2.0, 9.0])
def test_no_fit_based_failure_reason_is_ever_emitted(r2, n):
    result = classify(r2, n, True, DEFAULT_PARAMS)
    assert result.failure_reason is None
    r2_ok = math.isfinite(r2) and r2 >= DEFAULT_PARAMS.r2_threshold_good
    good = r2_ok and DEFAULT_PARAMS.n_min <= n <= DEFAULT_PARAMS.n_max
    assert result.affinity_class == ("Good" if good else "Partial")
