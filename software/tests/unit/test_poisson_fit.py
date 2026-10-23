"""Behavioral tests for poisson_fit — the Adams 2016 Appendix E estimator."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st
from scipy.optimize import minimize

from poisson_fit import baseline, exposure, library_fraction, predicted_mean_bin, profile_loglik

LOG_VALUES = 0.7 * np.arange(4)


def _brute_force_profile(R: np.ndarray, a: np.ndarray, k: float, ell: np.ndarray, m: float) -> float:
    """Oracle: max_q sum R log(k a q) - k a q, q on the simplex with sum q l = m, by SLSQP from many starts."""

    def nll(q: np.ndarray) -> float:
        x = k * a * np.maximum(q, 1e-300)
        return -(np.where(R > 0, R * np.log(x), 0.0) - x).sum()

    cons = [{"type": "eq", "fun": lambda q: q.sum() - 1}, {"type": "eq", "fun": lambda q: q @ ell - m}]
    best = np.inf
    rng = np.random.default_rng(0)
    for _ in range(12):
        q0 = rng.dirichlet(np.ones(len(R)))
        res = minimize(
            nll, q0, method="SLSQP", bounds=[(0, 1)] * len(R), constraints=cons, options={"ftol": 1e-14, "maxiter": 500}
        )
        # Tight: a constraint slip of 1e-7 can buy more likelihood than the 1e-7 relative tolerance below.
        feasible = abs(res.x.sum() - 1) < 1e-11 and abs(res.x @ ell - m) < 1e-11 and res.x.min() > -1e-12
        if feasible:
            best = min(best, res.fun)
    return -best


counts = st.lists(st.integers(0, 400), min_size=4, max_size=4).filter(lambda r: sum(r) > 0)


class TestProfileLoglik:
    # The KKT solve must never stop at a worse point than a general-purpose optimiser.
    @settings(max_examples=40, deadline=None)
    @given(
        R=counts,
        a=st.lists(st.floats(1e3, 1e6), min_size=4, max_size=4),
        m_frac=st.floats(0.05, 0.95),
    )
    def test_never_worse_than_a_general_optimiser(self, R, a, m_frac):
        R, a = np.array(R, float), np.array(a)
        k = (R / a).sum()  # the ratio estimate, as library_fraction gives for one concentration
        m = m_frac * LOG_VALUES[-1]
        ours = profile_loglik(R[None, None], a[None, None], np.array([k]), LOG_VALUES, np.array([m]))[0, 0, 0]
        oracle = _brute_force_profile(R, a, k, LOG_VALUES, m)
        assume(np.isfinite(oracle))  # the oracle found no feasible point: nothing to compare against
        assert ours >= oracle - 1e-7 * max(1.0, abs(oracle))

    # Exact value where it is known: at the mean the reads imply, with k = sum R/a, the optimum is
    # x = R, so the profile equals sum R log R - N. Above it would mean an infeasible q.
    @settings(max_examples=60, deadline=None)
    @given(R=counts, a=st.lists(st.floats(1e3, 1e6), min_size=4, max_size=4))
    def test_equals_the_closed_form_at_the_reads_implied_mean(self, R, a):
        R, a = np.array(R, float), np.array(a)
        q = (R / a) / (R / a).sum()
        m_hat = float(q @ LOG_VALUES)
        assume(0.01 < m_hat < LOG_VALUES[-1] - 0.01)
        ours = profile_loglik(R[None, None], a[None, None], np.array([(R / a).sum()]), LOG_VALUES, np.array([m_hat]))
        exact = (R[R > 0] * np.log(R[R > 0])).sum() - R.sum()
        assert ours[0, 0, 0] == pytest.approx(exact, abs=1e-7 * max(1.0, abs(exact)))

    # Without the curve constraint the best q is proportional to R/a; the profile over m must peak there.
    @settings(max_examples=40, deadline=None)
    @given(R=counts, a=st.lists(st.floats(1e3, 1e6), min_size=4, max_size=4))
    def test_peaks_at_the_reads_implied_mean(self, R, a):
        R, a = np.array(R, float), np.array(a)
        q = (R / a) / (R / a).sum()
        m_hat = float(q @ LOG_VALUES)
        grid = np.linspace(0.01, LOG_VALUES[-1] - 0.01, 200)
        ll = profile_loglik(R[None, None], a[None, None], np.array([(R / a).sum()]), LOG_VALUES, grid)[0, 0]
        assert (
            abs(grid[np.argmax(ll)] - m_hat) <= (grid[1] - grid[0]) + 1e-9 or not 0.01 < m_hat < LOG_VALUES[-1] - 0.01
        )

    def test_concentration_without_reads_carries_no_information(self):
        R = np.zeros((1, 1, 4))
        ll = profile_loglik(R, np.full((1, 1, 4), 1e4), np.array([1e-3]), LOG_VALUES, np.linspace(0.1, 2.0, 5))
        assert np.all(ll == 0.0)


class TestExposure:
    def test_bins_with_more_reads_than_cells_are_rescaled_to_cells(self):
        # Appendix E: each read must come from a different cell. 236,116 reads from 272 cells -> 272.
        reads = np.array([[1000.0, 236116.0]])
        cells = np.array([[5000.0, 272.0]])
        a, h = exposure(reads, cells)
        assert h[0, 1] == pytest.approx(272 / 236116)
        assert h[0, 0] == 1.0
        # a = rescaled reads * total cells / cells in bin
        assert a[0] == pytest.approx([1000 * 5272 / 5000, 272 * 5272 / 272])

    def test_sort_fractions_give_the_same_exposure_without_rescale(self):
        reads = np.array([[1000.0, 500.0]])
        cells = np.array([[3000.0, 1000.0]])
        a_abs, _ = exposure(reads, cells)
        a_frac, h = exposure(reads, cells / cells.sum(), absolute=False)
        assert np.all(h == 1.0)
        assert a_frac == pytest.approx(a_abs)

    def test_without_cells_every_bin_counts_equally(self):
        a, h = exposure(np.array([[10.0, 30.0]]), None)
        assert a[0] == pytest.approx([20.0, 60.0])
        assert np.all(h == 1.0)


def test_baseline_is_mean_linear_fluorescence_at_zero():
    # Two clonotypes: one all in bin 1 (value 1), one split evenly across bins 1 and 2 (1 and 10^0.7).
    reads0 = np.array([[100.0, 0, 0, 0], [50.0, 50.0, 0, 0]])
    b = baseline(reads0, np.ones(4), LOG_VALUES)
    assert b == pytest.approx((1.0 + (1.0 + 10**0.7) / 2) / 2)


def test_baseline_needs_reads_at_zero():
    assert baseline(np.zeros((2, 4)), np.ones(4), LOG_VALUES) is None


def test_library_fraction_is_geometric_mean_over_concentrations():
    reads = np.array([[[10.0, 0.0], [40.0, 0.0]]])
    a = np.ones((1, 2, 2)) * 10.0
    assert library_fraction(reads, a)[0] == pytest.approx(np.sqrt(1.0 * 4.0))


@settings(max_examples=50, deadline=None)
@given(log_kd=st.floats(-12, -2), amp=st.floats(0.1, 1e3), base=st.floats(1.0, 10.0))
def test_predicted_mean_bin_stays_within_the_bins(log_kd, amp, base):
    y = predicted_mean_bin(log_kd, amp, np.logspace(-10, -5, 11), base, 0.7, 4)
    assert np.all((y >= 1.0) & (y <= 4.0))
    assert np.all(np.diff(y) >= -1e-12)  # binding never decreases with concentration
