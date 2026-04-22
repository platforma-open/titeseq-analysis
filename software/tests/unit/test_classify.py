"""Behavioral tests for classify.py (R12 truth table + boundaries)."""

from __future__ import annotations

import pytest

from classify import classify
from constants import DEFAULT_PARAMS, FitParams


# Truth table: covers every spec row + all boundaries + negative R² passthrough.
@pytest.mark.parametrize(
    "r2, n, converged, expected_class, expected_reason",
    [
        # Row 1: R² >= Good, n in range → Good
        (0.90, 1.0, True, "Good", None),
        (0.80, 1.0, True, "Good", None),  # boundary: R² == r2_threshold_good (inclusive)
        # Row 2: R² >= Good, n out of range → Partial (downgrade, not Failed)
        (0.95, 0.4, True, "Partial", None),
        (0.95, 2.1, True, "Partial", None),
        (0.95, 0.5, True, "Good", None),  # boundary: n == n_min (inclusive)
        (0.95, 2.0, True, "Good", None),  # boundary: n == n_max (inclusive)
        # Row 3: Failed <= R² < Good, n in range → Partial
        (0.50, 1.0, True, "Partial", None),  # boundary: R² == r2_threshold_failed
        (0.79, 1.0, True, "Partial", None),
        # Row 4: Failed <= R² < Good, n out of range → Failed (n_out_of_range)
        (0.70, 0.3, True, "Failed", "n_out_of_range"),
        (0.70, 2.5, True, "Failed", "n_out_of_range"),
        # Row 5: R² < r2_threshold_failed → Failed (low_r2)
        (0.49, 1.0, True, "Failed", "low_r2"),
        (-0.5, 1.0, True, "Failed", "low_r2"),  # negative R² path — classification still fires
        # Did not converge → Failed (convergence_failure), regardless of R²/n
        (0.0, 0.0, False, "Failed", "convergence_failure"),
        (None, None, False, "Failed", "convergence_failure"),
    ],
)
def test_classification_truth_table(r2, n, converged, expected_class, expected_reason):
    result = classify(r2, n, converged, DEFAULT_PARAMS)
    assert result.affinity_class == expected_class
    assert result.failure_reason == expected_reason


# Non-default thresholds (e.g. multimeric-antigen scenario) exercise the same truth table.
# Guards against hardcoded constants leaking into classify().
def test_classification_with_non_default_thresholds():
    strict = FitParams(r2_threshold_good=0.95, r2_threshold_failed=0.7, n_min=0.8, n_max=3.0)
    # R² = 0.93, n=1.5: in Partial zone under strict thresholds (would be Good with defaults).
    result = classify(0.93, 1.5, True, strict)
    assert result.affinity_class == "Partial"
    assert result.failure_reason is None

    # R² = 0.96, n=0.5: n_min is 0.8 so n is OOR → Partial (downgrade).
    result = classify(0.96, 0.5, True, strict)
    assert result.affinity_class == "Partial"


# Precomputed upstream gates (insufficient_*, non_monotonic_signal) win over fit-based reasons.
@pytest.mark.parametrize(
    "precomputed",
    [
        "insufficient_reads",
        "insufficient_points",
        "non_monotonic_signal",
    ],
)
def test_precomputed_reason_wins(precomputed):
    result = classify(0.99, 1.0, True, DEFAULT_PARAMS, precomputed_reason=precomputed)
    assert result.affinity_class == "Failed"
    assert result.failure_reason == precomputed
