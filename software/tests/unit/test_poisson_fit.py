"""Behavioral tests for poisson_fit — the Adams 2016 Appendix E estimator."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import assume, example, given, settings
from hypothesis import strategies as st
from scipy.optimize import minimize

from poisson_fit import (
    BASELINE_CONTAMINATION_RATIO,
    BASELINE_MEDIAN_DIVISOR,
    BASELINE_MIN_READS,
    baseline,
    exposure,
    fit_titration,
    library_fraction,
    predicted_mean_bin,
    predicted_mean_bin_measured,
    profile_loglik,
    read_weighted_median,
)

LOG_VALUES = 0.7 * np.arange(4)
# profile_loglik and fit_titration take one row of bin log-fluorescence per concentration.
ONE_CONC = np.broadcast_to(LOG_VALUES, (1, LOG_VALUES.size))
BIN_LABELS = np.arange(1, 5)


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
        ours = profile_loglik(R[None, None], a[None, None], np.array([k]), ONE_CONC, np.array([m]))[0, 0, 0]
        oracle = _brute_force_profile(R, a, k, LOG_VALUES, m)
        assume(np.isfinite(oracle))  # the oracle found no feasible point: nothing to compare against
        assert ours >= oracle - 1e-7 * max(1.0, abs(oracle))

    # Bins that share one exposure are common in real data (equal cells sorted per bin). When an empty bin
    # shares it, the one-empty-bin root sits on a pole at nu = 0, and the bisection approaches it through
    # subnormal floats without stopping early. The pinned example is such a case.
    @settings(max_examples=40, deadline=None)
    @example(R=[300, 70, 0, 50], exposures=[None, 130.0, None, None], shared=1000.0, m_frac=1.0 / 2.1)
    @given(
        R=counts.filter(lambda r: 0 in r),
        exposures=st.lists(st.none() | st.floats(1e3, 1e6), min_size=4, max_size=4),
        shared=st.floats(1e3, 1e6),
        m_frac=st.floats(0.05, 0.95),
    )
    def test_shared_exposure_never_worse_than_a_general_optimiser(self, R, exposures, shared, m_frac):
        # Every empty bin, and every bin drawn as None, gets the shared exposure.
        a = np.array([shared if r == 0 or e is None else e for r, e in zip(R, exposures)])
        R = np.array(R, float)
        k = (R / a).sum()
        m = m_frac * LOG_VALUES[-1]
        ours = profile_loglik(R[None, None], a[None, None], np.array([k]), ONE_CONC, np.array([m]))[0, 0, 0]
        oracle = _brute_force_profile(R, a, k, LOG_VALUES, m)
        assume(np.isfinite(oracle))
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
        ours = profile_loglik(R[None, None], a[None, None], np.array([(R / a).sum()]), ONE_CONC, np.array([m_hat]))
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
        ll = profile_loglik(R[None, None], a[None, None], np.array([(R / a).sum()]), ONE_CONC, grid)[0, 0]
        assert (
            abs(grid[np.argmax(ll)] - m_hat) <= (grid[1] - grid[0]) + 1e-9 or not 0.01 < m_hat < LOG_VALUES[-1] - 0.01
        )

    def test_concentration_without_reads_carries_no_information(self):
        R = np.zeros((1, 1, 4))
        ll = profile_loglik(R, np.full((1, 1, 4), 1e4), np.array([1e-3]), ONE_CONC, np.linspace(0.1, 2.0, 5))
        assert np.all(ll == 0.0)

    # Each concentration is profiled against its own row of bin fluorescence: identical reads under
    # two different gate scales give two different likelihoods.
    def test_each_concentration_uses_its_own_bin_fluorescence(self):
        reads = np.array([[[100.0, 50.0, 10.0, 1.0], [100.0, 50.0, 10.0, 1.0]]])
        a = np.full((1, 2, 4), 1e4)
        k = np.array([1e-3])
        m_grid = np.array([0.5, 0.8])
        narrow = 0.4 * np.arange(4)
        both = profile_loglik(reads, a, k, np.vstack([LOG_VALUES, narrow]), m_grid)

        first = profile_loglik(reads[:, :1], a[:, :1], k, ONE_CONC, m_grid)
        second = profile_loglik(reads[:, 1:], a[:, 1:], k, narrow[None, :], m_grid)
        assert both[:, 0] == pytest.approx(first[:, 0])
        assert both[:, 1] == pytest.approx(second[:, 0])
        assert not np.allclose(both[:, 0], both[:, 1])


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


class TestReadWeightedMedian:
    @pytest.mark.parametrize(
        "sizes, expected",
        [
            ([], 0.0),
            ([0.0, 0.0, 0.0], 0.0),
            ([5.0], 5.0),
            ([1.0] * 10, 1.0),
            ([10.0, 20.0, 30.0, 40.0], 30.0),  # half of 100 reads lie below 30
            # The plain median of these is 1; the reads live in the one deep clonotype.
            ([1.0, 1.0, 1.0, 1.0, 1000.0], 1000.0),
        ],
    )
    def test_pinned(self, sizes, expected):
        assert read_weighted_median(np.array(sizes, dtype=float)) == pytest.approx(expected)

    def test_is_invariant_to_order(self):
        sizes = np.array([7.0, 1.0, 900.0, 3.0, 40.0])
        assert read_weighted_median(sizes) == read_weighted_median(sizes[::-1])

    # A crowd of one-read clonotypes carrying less read mass than the real library collapses the
    # plain median to 1 and leaves this one where it was.
    def test_unaffected_by_a_crowd_of_one_read_clonotypes(self):
        real = np.array([800.0, 1000.0, 1200.0])
        with_errors = np.concatenate([real, np.ones(1_000)])
        assert np.median(with_errors) == 1.0
        assert read_weighted_median(with_errors) == read_weighted_median(real) == 1000.0


class TestBaseline:
    def test_is_mean_linear_fluorescence_at_zero(self):
        # Two clonotypes: one all in bin 1 (value 1), one split evenly across bins 1 and 2 (1 and 10^0.7).
        reads0 = np.array([[100.0, 0, 0, 0], [50.0, 50.0, 0, 0]])
        est = baseline(reads0, np.ones(4), LOG_VALUES, reads0)
        assert est.value == pytest.approx((1.0 + (1.0 + 10**0.7) / 2) / 2)
        assert (est.n_kept, est.n_total) == (2, 2)
        assert est.read_floor == BASELINE_MIN_READS

    def test_value_is_none_only_without_any_signal_at_zero(self):
        est = baseline(np.zeros((2, 4)), np.ones(4), LOG_VALUES, np.zeros((2, 4)))
        assert est.value is None
        assert (est.n_kept, est.n_total) == (0, 0)
        assert est.contamination_ratio is None

    # The floor keeps shallow clonotypes out of the mean: each would otherwise carry equal weight
    # whatever its depth, and a one-read clonotype in the top bin implies a fluorescence near 10^2.1.
    def test_shallow_clonotypes_are_excluded(self):
        deep = np.array([[1000.0, 0.0, 0.0, 0.0]])
        shallow = np.tile(np.array([[0.0, 0.0, 0.0, 1.0]]), (20, 1))
        reads0 = np.vstack([deep, shallow])
        est = baseline(reads0, np.ones(4), LOG_VALUES, reads0)
        assert est.read_floor == BASELINE_MIN_READS
        assert (est.n_kept, est.n_total) == (1, 21)
        assert est.value == pytest.approx(1.0)

    def test_floor_scales_with_the_read_weighted_median(self):
        # Median clonotype size 6000 → floor 200, above BASELINE_MIN_READS.
        reads0 = np.array([[6000.0, 0, 0, 0], [150.0, 0, 0, 0], [6000.0, 0, 0, 0]])
        est = baseline(reads0, np.ones(4), LOG_VALUES, reads0)
        assert est.read_floor == pytest.approx(6000.0 / BASELINE_MEDIAN_DIVISOR)
        assert est.read_floor > BASELINE_MIN_READS
        assert (est.n_kept, est.n_total) == (2, 3)

    # A library too shallow for any clonotype to clear the floor still needs an answer.
    def test_floor_drops_to_zero_when_nothing_clears_it(self):
        reads0 = np.tile(np.array([[1.0, 0.0, 0.0, 0.0]]), (20, 1))
        est = baseline(reads0, np.ones(4), LOG_VALUES, reads0)
        assert est.read_floor == 0.0
        assert (est.n_kept, est.n_total) == (20, 20)
        assert est.value == pytest.approx(1.0)

    # The rescale shrinks the high bins by two to three orders of magnitude. Flooring on the
    # rescaled reads would delete exactly the high-fluorescence clonotypes the estimate needs.
    def test_floor_is_applied_to_raw_reads_not_rescaled_ones(self):
        raw = np.array([[0.0, 0.0, 0.0, 5000.0], [200.0, 0.0, 0.0, 0.0]])
        rescaled = np.array([[0.0, 0.0, 0.0, 5.0], [200.0, 0.0, 0.0, 0.0]])
        est = baseline(rescaled, np.ones(4), LOG_VALUES, raw)
        assert est.n_kept == 2
        assert est.value == pytest.approx((10**2.1 + 1.0) / 2)
        # The high-fluorescence clonotype keeps 5 reads after the rescale: a floor on those would
        # drop it and leave a value of 1.0.
        assert rescaled.sum(1)[0] < BASELINE_MIN_READS <= est.read_floor

    def test_contamination_ratio_is_one_on_a_clean_zero_molar_set(self):
        reads0 = np.array([[500.0, 0, 0, 0], [250.0, 250.0, 0, 0]])
        est = baseline(reads0, np.ones(4), LOG_VALUES, reads0)
        assert est.contamination_ratio == pytest.approx(1.0)
        assert est.contamination_ratio < BASELINE_CONTAMINATION_RATIO

    # Shallow clonotypes in the high bins pull the unfloored clonotype mean far above the pooled one.
    def test_contamination_ratio_flags_a_crowd_of_shallow_high_bin_clonotypes(self):
        deep = np.array([[1000.0, 0.0, 0.0, 0.0]])
        shallow = np.tile(np.array([[0.0, 0.0, 0.0, 1.0]]), (20, 1))
        reads0 = np.vstack([deep, shallow])
        est = baseline(reads0, np.ones(4), LOG_VALUES, reads0)
        assert est.contamination_ratio > BASELINE_CONTAMINATION_RATIO


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


class TestPredictedMeanBinMeasured:
    CONCS = np.logspace(-10, -5, 11)
    SCALE = np.broadcast_to(LOG_VALUES, (CONCS.size, LOG_VALUES.size))

    # Against the uniform scale the fixed-scale and measured-scale predictions are the same curve.
    @settings(max_examples=50, deadline=None)
    @given(log_kd=st.floats(-12, -2), amp=st.floats(0.1, 1e3), base=st.floats(1.0, 10.0))
    def test_agrees_with_the_fixed_scale_when_the_scale_is_uniform(self, log_kd, amp, base):
        fixed = predicted_mean_bin(log_kd, amp, self.CONCS, base, 0.7, 4)
        measured = predicted_mean_bin_measured(log_kd, amp, self.CONCS, base, self.SCALE, BIN_LABELS)
        assert measured == pytest.approx(fixed, abs=1e-9)

    # A gate that sits higher at high concentration puts the same curve in a lower bin there.
    def test_follows_each_concentrations_own_gates(self):
        drifting = np.vstack([LOG_VALUES + 0.5 * i / (self.CONCS.size - 1) for i in range(self.CONCS.size)])
        uniform = predicted_mean_bin_measured(-8.0, 50.0, self.CONCS, 2.5, self.SCALE, BIN_LABELS)
        drifted = predicted_mean_bin_measured(-8.0, 50.0, self.CONCS, 2.5, drifting, BIN_LABELS)
        assert drifted[-1] < uniform[-1]
        assert drifted[0] == pytest.approx(uniform[0])

    def test_stays_within_the_bin_labels(self):
        y = predicted_mean_bin_measured(-12.0, 1e4, self.CONCS, 2.5, self.SCALE, BIN_LABELS)
        assert np.all((y >= 1.0) & (y <= 4.0))


def test_fit_titration_recovers_a_planted_kd():
    concentrations = np.array([0.0] + [10 ** (-10 + 0.5 * i) for i in range(10)])
    baseline_fluorescence, amplitude, true_log10_kd = 2.5, 50.0, -8.0
    log_values = np.broadcast_to(LOG_VALUES, (concentrations.size, LOG_VALUES.size))
    # Noiseless reads: bin probabilities are a Gaussian (sd 0.3 decades) about the model mean.
    m = np.log10(baseline_fluorescence + amplitude * concentrations / (concentrations + 10.0**true_log10_kd))
    q = np.exp(-0.5 * ((LOG_VALUES[None, :] - m[:, None]) / 0.3) ** 2)
    q /= q.sum(1, keepdims=True)
    reads = (q * 2000.0)[None, :, :]
    a = np.full((1, concentrations.size, 4), 1e5)

    fit = fit_titration(reads, a, concentrations, log_values, baseline_fluorescence, (1.0, 1000.0))
    assert fit.log10_kd[0] == pytest.approx(true_log10_kd, abs=0.1)
    assert fit.amplitude[0] == pytest.approx(amplitude, rel=0.2)
    assert fit.log10_kd_sd[0] < 0.5
