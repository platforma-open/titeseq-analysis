"""Weighted R² — and a guard against the two copies of it drifting apart.

`hill_fit` carries its own copy so the v2.4.0 path is untouched by this change. The equality test
below is what stops the two from diverging unnoticed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from fit_stats import weighted_r2
from hill_fit import weighted_r2 as hill_weighted_r2

Y = np.array([1.0, 2.0, 3.0, 4.0, 5.0])


def test_a_perfect_fit_scores_one():
    assert weighted_r2(Y, Y, np.ones(5)) == pytest.approx(1.0)


def test_unit_weights_give_the_ordinary_r2():
    y_hat = np.array([1.1, 2.1, 2.9, 4.2, 4.8])
    expected = 1.0 - ((Y - y_hat) ** 2).sum() / ((Y - Y.mean()) ** 2).sum()
    assert weighted_r2(Y, y_hat, np.ones(5)) == pytest.approx(expected)


def test_a_worse_than_mean_fit_is_not_clamped():
    assert weighted_r2(Y, np.full(5, 100.0), np.ones(5)) < 0.0


def test_zero_weight_points_contribute_to_neither_sum():
    w = np.array([1.0, 1.0, 1.0, 1.0, 0.0])
    y_hat = np.array([1.0, 2.0, 3.0, 4.0, 99.0])
    assert weighted_r2(Y, y_hat, w) == pytest.approx(1.0)


def test_no_weight_and_no_spread_both_give_nan():
    assert math.isnan(weighted_r2(Y, Y, np.zeros(5)))
    assert math.isnan(weighted_r2(np.full(5, 2.0), np.full(5, 2.0), np.ones(5)))


@pytest.mark.parametrize(
    "y_hat,w",
    [
        (np.array([1.1, 2.1, 2.9, 4.2, 4.8]), np.ones(5)),
        (np.array([0.5, 2.5, 3.5, 3.5, 5.5]), np.array([10.0, 1.0, 3.0, 0.0, 7.0])),
        (np.full(5, 100.0), np.array([1.0, 2.0, 3.0, 4.0, 5.0])),
    ],
)
def test_the_two_copies_agree(y_hat, w):
    assert weighted_r2(Y, y_hat, w) == pytest.approx(hill_weighted_r2(Y, y_hat, w))
