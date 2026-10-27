"""K_D inference from read counts — Adams, Mora, Walczak, Kinney 2016 (eLife 5:e23156), Appendix E.

Model, per clonotype s and concentration c:
    mean log10 fluorescence  m_sc = log10(B + A_s * c / (c + K_s))          (E1, E2; Hill n = 1)
    expected reads           r_bsc = a_bc * k_s * q_b|sc,  a_bc = T_bc * C_c / C_bc   (E4)
    with q_b|sc >= 0, sum_b q_b|sc = 1, sum_b q_b|sc * l_b = m_sc                  (E5)
where l_bc is the bin's log10 fluorescence, T_bc the reads in sample (b, c), C_bc the cells sorted
into it, and k_s the clonotype's library fraction. Absent measured values l_bc = bin_width *
(bin_label - 1), one fixed scale for every concentration; with them it is the sort report's own mean
per gate, which moves with concentration because at low concentration the cells in the high gates
sit at the bottom of the gate.
Reads are Poisson (E3). The bin probabilities q are profiled out; K and A are found on a grid.

k_s is the ratio estimate (geometric mean over concentrations of sum_b R_b / a_b), shared by all
concentrations. For fixed k and m the profile over q is concave with two linear constraints; it is
solved exactly from its KKT conditions, vectorised over clonotypes, concentrations and a grid of m.
The likelihood at the m a (K, A) pair predicts is interpolated from that grid.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_NEWTON_STEPS = 60
_BISECTION_STEPS = 200
# Feasibility tolerance on the profile's constraints (probabilities sum to 1, mean equals m).
_TOL = 1e-7
_BACKTRACK_STEPS = 40
# Spacing of the K grid in decades. Adams: 201 points over [1e-10, 1e-4].
KD_GRID_STEP = 0.03
AMPLITUDE_GRID_POINTS = 161
M_GRID_POINTS = 120

# Read floor for a clonotype to count toward the 0 M baseline: max(100, M/30), M the read-weighted
# median clonotype size. The unfloored mean converges to a computable wrong number as
# sequencing-error clonotypes dominate, because each carries equal weight whatever its depth. A
# one-read clonotype implies a fluorescence near 10 against a correct 2.59.
BASELINE_MIN_READS = 100.0
BASELINE_MEDIAN_DIVISOR = 30.0
# Unfloored clonotype mean over the pooled baseline, above which the 0 M set looks contaminated.
# Clean Adams 2016 replicates span 0.79 to 1.26; one-fold contamination starts at 2.03.
BASELINE_CONTAMINATION_RATIO = 1.5

# Largest spread, in decades, the measured top gate may show across concentrations before the
# measured scale is refused. The model gives one fluorescence per (gate, concentration), shared by
# every clonotype. That holds while a gate's occupants are representative. Once the top gate's own
# mean collapses at low concentration — because the cells sitting there are background rather than
# bound — a saturated binder's mean lies above the top gate, the profile likelihood reads that as
# impossible, and no tight binder can be fitted at all.
BIN_FLUORESCENCE_MAX_TOP_SPREAD = 0.5


@dataclass(frozen=True)
class BaselineEstimate:
    """B plus what it was computed from, so the pipeline can report and warn on it."""

    value: float | None
    n_kept: int
    n_total: int
    read_floor: float
    # Unfloored clonotype mean divided by the pooled baseline; None when either is unavailable.
    contamination_ratio: float | None


@dataclass(frozen=True)
class GridFit:
    log10_kd: np.ndarray  # (V,) best K on the grid, unclipped
    amplitude: np.ndarray  # (V,) best A, in units of bin-1 fluorescence
    log10_kd_sd: np.ndarray  # (V,) posterior RMS distance of log10 K from the best fit, flat prior on the grid
    log10_kd_grid: np.ndarray


def exposure(
    sample_reads: np.ndarray, sample_cells: np.ndarray | None, *, absolute: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    """a_bc = T_bc * C_c / C_bc, plus the read rescale h_bc = min(1, C_bc / T_bc) (Appendix E, after E4).

    sample_reads, sample_cells: (C, B) totals per sample. `absolute=False` means sample_cells are
    sort fractions: the ratio C_c / C_bc is still right, but reads cannot be compared with cells, so no
    rescale. With no cell counts every bin is treated as holding the same number of cells.
    """
    T = np.asarray(sample_reads, dtype=float)
    if sample_cells is None:
        return T * T.shape[1], np.ones_like(T)
    C = np.asarray(sample_cells, dtype=float)
    h = np.where((T > C) & (C > 0) & absolute, C / np.where(T > 0, T, 1.0), 1.0)
    Th = T * h
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(C > 0, Th * C.sum(1, keepdims=True) / C, 0.0)
    return a, h


def library_fraction(reads: np.ndarray, a: np.ndarray) -> np.ndarray:
    """k_s: geometric mean over concentrations of sum_b R_bsc / a_bc (Adams ratio_x)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        kc = np.where(a > 0, reads / np.where(a > 0, a, 1.0), 0.0).sum(-1)
        lk = np.where(kc > 0, np.log(np.where(kc > 0, kc, 1.0)), np.nan)
    cnt = np.isfinite(lk).sum(1)
    mean = np.where(cnt > 0, np.nansum(np.where(np.isfinite(lk), lk, 0.0), 1) / np.maximum(cnt, 1), np.log(1e-12))
    return np.exp(mean)


def _implied_fluorescence(reads0: np.ndarray, a0: np.ndarray, log_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(V,) linear fluorescence each clonotype's 0 M bin probabilities imply, plus which are defined."""
    with np.errstate(divide="ignore", invalid="ignore"):
        q = np.where(a0 > 0, reads0 / np.where(a0 > 0, a0, 1.0), 0.0)
    tot = q.sum(1)
    defined = tot > 0
    lin = np.full(reads0.shape[0], np.nan)
    if defined.any():
        lin[defined] = ((q[defined] / tot[defined, None]) * 10.0**log_values).sum(1)
    return lin, defined


def read_weighted_median(sizes: np.ndarray) -> float:
    """Clonotype size at which half the reads lie below it — the N50 of clonotype read counts.

    The plain median tracks how many error clonotypes there are; this tracks how deeply the real
    library was sequenced, which is what the floor has to scale with.
    """
    s = np.sort(np.asarray(sizes, dtype=float))
    s = s[s > 0]
    if s.size == 0:
        return 0.0
    c = np.cumsum(s)
    return float(s[int(np.searchsorted(c, c[-1] / 2.0))])


def baseline(reads0: np.ndarray, a0: np.ndarray, log_values: np.ndarray, raw_reads0: np.ndarray) -> BaselineEstimate:
    """B: mean implied 0 M fluorescence over clonotypes deep enough to carry the estimate.

    reads0: (V, B) exposure-rescaled reads at 0 M; a0: (B,) exposure at 0 M; raw_reads0: (V, B) the
    same reads before the rescale. Appendix E sets B to the mean fluorescence of cells labelled at
    0 M; in bin-1 units this reproduces the authors' values to 3 decimals.

    The floor is applied to the RAW reads, never the rescaled ones. The rescale shrinks the high
    bins by two to three orders of magnitude, so flooring on it removes exactly the
    high-fluorescence clonotypes the estimate needs.
    """
    lin, defined = _implied_fluorescence(reads0, a0, log_values)
    raw_total = np.asarray(raw_reads0, dtype=float).sum(1)
    read_floor = max(BASELINE_MIN_READS, read_weighted_median(raw_total) / BASELINE_MEDIAN_DIVISOR)

    keep = defined & (raw_total >= read_floor)
    # A library sequenced too shallowly for any clonotype to clear the floor still needs an answer.
    if not keep.any():
        keep, read_floor = defined, 0.0
    if not keep.any():
        return BaselineEstimate(None, 0, 0, read_floor, None)

    pooled_lin, pooled_ok = _implied_fluorescence(reads0.sum(0, keepdims=True), a0, log_values)
    unfloored = float(np.nanmean(lin[defined]))
    ratio = unfloored / float(pooled_lin[0]) if bool(pooled_ok[0]) and float(pooled_lin[0]) > 0 else None
    return BaselineEstimate(float(np.nanmean(lin[keep])), int(keep.sum()), int(defined.sum()), read_floor, ratio)


def _loglik(R: np.ndarray, x: np.ndarray) -> np.ndarray:
    """sum_b R_b log x_b - x_b, with expected reads x (Poisson, without the log R! constant)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return (np.where(R > 0, R * np.log(np.where(x > 0, x, 1.0)), 0.0) - x).sum(-1)


def _interior(
    R: np.ndarray, ka: np.ndarray, L: np.ndarray, m: np.ndarray, filled: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """No empty bin carries mass: damped Newton on the dual over the bins with reads.

    Only points whose m lies strictly inside the range of the filled bins can be interior; the rest
    are skipped. Converged points leave the working set, so the cost follows the hard points only.
    """
    shape, B = m.shape, R.shape[-1]
    Rf, kaf, Lf, ff = R.reshape(-1, B), ka.reshape(-1, B), L.reshape(-1, B), filled.reshape(-1, B)
    mf = m.reshape(-1)
    lo = np.where(ff, Lf, np.inf).min(-1)
    hi = np.where(ff, Lf, -np.inf).max(-1)
    lam = Rf.sum(-1) + 1.0
    nu = np.zeros_like(mf)
    act = np.flatnonzero((lo < mf) & (mf < hi))

    def dual(i: np.ndarray, lam_: np.ndarray, nu_: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        s_ = kaf[i] + lam_[:, None] + nu_[:, None] * Lf[i]
        bad = (ff[i] & (s_ <= 0)).any(-1)
        val = -np.where(ff[i], Rf[i] * np.log(np.where(s_ > 0, s_, 1.0)), 0.0).sum(-1) + lam_ + nu_ * mf[i]
        return np.where(bad, np.inf, val), s_

    for _ in range(_NEWTON_STEPS):
        if act.size == 0:
            break
        la, na = lam[act], nu[act]
        D, s = dual(act, la, na)
        w = np.where(ff[act], Rf[act] / s, 0.0)
        g1 = 1.0 - w.sum(-1)
        g2 = mf[act] - (w * Lf[act]).sum(-1)
        conv = (np.abs(g1) < 1e-13) & (np.abs(g2) < 1e-13)
        w2 = np.where(ff[act], w / s, 0.0)
        h11, h12, h22 = w2.sum(-1), (w2 * Lf[act]).sum(-1), (w2 * Lf[act] ** 2).sum(-1)
        det = h11 * h22 - h12**2
        det = np.where(det > 0, det, np.inf)
        d1 = -(h22 * g1 - h12 * g2) / det
        d2 = -(-h12 * g1 + h11 * g2) / det
        t = np.ones(act.size)
        done = conv.copy()
        moved = np.zeros(act.size, dtype=bool)
        for _ in range(_BACKTRACK_STEPS):
            cl, cn = la + t * d1, na + t * d2
            cD, _ = dual(act, cl, cn)
            ok = (~done) & (cD <= D + 1e-12 * np.abs(D))  # near the optimum D is flat to rounding
            la, na = np.where(ok, cl, la), np.where(ok, cn, na)
            done |= ok
            moved |= ok
            if done.all():
                break
            t = np.where(done, t, 0.5 * t)
        lam[act], nu[act] = la, na
        act = act[~conv & moved]

    s = kaf + lam[:, None] + nu[:, None] * Lf
    q = np.where(ff, Rf / np.where(s > 0, s, 1.0), 0.0)
    feasible = (lo < mf) & (mf < hi)
    feasible &= (np.abs(q.sum(-1) - 1.0) < _TOL) & (np.abs((q * Lf).sum(-1) - mf) < _TOL)
    feasible &= np.where(ff, s > 0, s >= -_TOL * np.abs(kaf)).all(-1)
    return q.reshape(R.shape), feasible.reshape(shape)


def _one_empty_bin(
    R: np.ndarray, ka: np.ndarray, L: np.ndarray, m: np.ndarray, filled: np.ndarray, z: int
) -> tuple[np.ndarray, np.ndarray]:
    """Empty bin z carries mass (its dual slack is 0): a monotone 1-D root in nu."""
    c = ka - ka[..., z : z + 1]
    d = L - L[..., z : z + 1]
    rhs = m - L[..., z]
    with np.errstate(divide="ignore", invalid="ignore"):
        bound = -c / d
    lo = np.where(filled & (d > 0), bound, -np.inf).max(-1)
    hi = np.where(filled & (d < 0), bound, np.inf).min(-1)
    # Each term R_b d_b / (c_b + nu d_b) is at most R_b / (nu - lo) above the pole lo, and at least
    # -R_b / (hi - nu) below the pole hi, so the root lies within sum(R) / |rhs| of the one finite
    # end. Twice that leaves room for rounding.
    with np.errstate(divide="ignore", invalid="ignore"):
        reach = 2.0 * R.sum(-1) / np.abs(rhs)
    one_sided_hi = np.isfinite(lo) & ~np.isfinite(hi) & (rhs > 0)
    one_sided_lo = np.isfinite(hi) & ~np.isfinite(lo) & (rhs < 0)
    lo, hi = np.where(one_sided_lo, hi - reach, lo), np.where(one_sided_hi, lo + reach, hi)
    span = 1e15 * (1.0 + np.abs(ka).max(-1) + R.sum(-1))
    lo = np.where(np.isfinite(lo), lo, np.minimum(hi, 0.0) - span)
    hi = np.where(np.isfinite(hi), hi, np.maximum(lo, 0.0) + span)

    def g(nu: np.ndarray) -> np.ndarray:
        den = c + nu[..., None] * d
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(filled, R * d / np.where(filled, den, 1.0), 0.0).sum(-1) - rhs

    for _ in range(_BISECTION_STEPS):  # g decreases in nu
        mid = 0.5 * (lo + hi)
        up = g(mid) > 0
        new_lo, new_hi = np.where(up, mid, lo), np.where(up, hi, mid)
        # Unchanged brackets repeat on every later step: stop once none moves.
        if not ((new_lo != lo) | (new_hi != hi)).any():
            break
        lo, hi = new_lo, new_hi
    nu = 0.5 * (lo + hi)
    slack = c + nu[..., None] * d
    q = np.where(filled, R / np.where(filled & (slack > 0), slack, np.inf), 0.0)
    qz = 1.0 - q.sum(-1)
    q[..., z] = qz
    feasible = (qz >= -_TOL) & (np.abs((q * L).sum(-1) - m) < _TOL)
    feasible &= np.where(filled, slack > 0, slack >= -_TOL * np.abs(ka)).all(-1)
    return np.maximum(q, 0.0), feasible


def _two_empty_bins(
    R: np.ndarray, ka: np.ndarray, L: np.ndarray, m: np.ndarray, filled: np.ndarray, z1: int, z2: int
) -> tuple[np.ndarray, np.ndarray]:
    """Empty bins z1, z2 carry mass: both slacks are 0, which fixes the duals; then a 2x2 solve."""
    nu = -(ka[..., z1] - ka[..., z2]) / (L[..., z1] - L[..., z2])
    lam = -ka[..., z1] - nu * L[..., z1]
    slack = ka + lam[..., None] + nu[..., None] * L
    q = np.where(filled, R / np.where(filled & (slack > 0), slack, np.inf), 0.0)
    r0 = 1.0 - q.sum(-1)
    r1 = m - (q * L).sum(-1)
    q2 = (r1 - L[..., z1] * r0) / (L[..., z2] - L[..., z1])
    q[..., z2] = q2
    q[..., z1] = r0 - q2
    feasible = (q[..., z1] >= -_TOL) & (q2 >= -_TOL)
    others = np.ones(L.shape[-1], dtype=bool)
    others[[z1, z2]] = False
    feasible &= np.where(filled, slack > 0, np.where(others, slack >= -_TOL * np.abs(ka), True)).all(-1)
    return np.maximum(q, 0.0), feasible


def profile_loglik(
    reads: np.ndarray, a: np.ndarray, k: np.ndarray, log_values: np.ndarray, m_grid: np.ndarray
) -> np.ndarray:
    """max over q of sum_b R_b log(k a_b q_b) - k a_b q_b, for every m in m_grid. Returns (V, C, M).

    reads, a: (V, C, B); k: (V,); log_values: (C, B). Rows with no reads give 0 (they carry no
    information on m).
    Solved exactly through the KKT conditions: either no empty bin carries probability (Newton on the
    two-variable dual), or one or two empty bins do (their dual slack is zero, which leaves a 1-D root
    or a closed form). Every feasible candidate is optimal; the best is kept.
    """
    shape = reads.shape[:2] + (m_grid.size,)
    B = reads.shape[-1]
    R = np.broadcast_to(reads[:, :, None, :], shape + (B,))
    ka = np.broadcast_to((k[:, None, None] * a)[:, :, None, :], shape + (B,))
    L = np.broadcast_to(np.asarray(log_values, dtype=float)[None, :, None, :], shape + (B,))
    m = np.broadcast_to(m_grid[None, None, :], shape)
    filled = R > 0
    q, ok = _interior(R, ka, L, m, filled)
    flat = np.where(ok, _loglik(R, ka * q), -np.inf).reshape(-1)
    # Empty-bin cases only where those bins are empty: solve on the flattened subset, scatter back.
    Rf, kaf, Lf, ff, mf = R.reshape(-1, B), ka.reshape(-1, B), L.reshape(-1, B), filled.reshape(-1, B), m.reshape(-1)
    for z in range(B):
        idx = np.flatnonzero(~ff[:, z])
        if idx.size:
            q, ok = _one_empty_bin(Rf[idx], kaf[idx], Lf[idx], mf[idx], ff[idx], z)
            flat[idx] = np.where(ok, np.maximum(flat[idx], _loglik(Rf[idx], kaf[idx] * q)), flat[idx])
    for z1 in range(B):
        for z2 in range(z1 + 1, B):
            idx = np.flatnonzero(~ff[:, z1] & ~ff[:, z2])
            if idx.size:
                q, ok = _two_empty_bins(Rf[idx], kaf[idx], Lf[idx], mf[idx], ff[idx], z1, z2)
                flat[idx] = np.where(ok, np.maximum(flat[idx], _loglik(Rf[idx], kaf[idx] * q)), flat[idx])
    return np.where(reads.sum(-1)[:, :, None] > 0, flat.reshape(shape), 0.0)


def fit_titration(
    reads: np.ndarray,
    a: np.ndarray,
    concentrations: np.ndarray,
    log_values: np.ndarray,
    baseline_fluorescence: float,
    amplitude_bounds: tuple[float, float],
    *,
    chunk: int = 25,
) -> GridFit:
    """Grid maximum likelihood of (K, A) for every clonotype.

    reads, a: (V, C, B) for all concentrations including 0 M (the 0 M row predicts m = log10 B).
    concentrations: (C,); log_values: (C, B). A model mean outside the bins' own fluorescence range
    contributes nothing to the likelihood, as in the authors' code, where it falls outside the
    interpolation grid.
    """
    conc = np.asarray(concentrations, dtype=float)
    log_values = np.asarray(log_values, dtype=float)
    pos = conc[conc > 0]
    lo, hi = np.log10(pos.min()) - 0.5, np.log10(pos.max()) + 1.0
    logk = np.linspace(lo, hi, int(round((hi - lo) / KD_GRID_STEP)) + 1)
    A = np.logspace(np.log10(amplitude_bounds[0]), np.log10(amplitude_bounds[1]), AMPLITUDE_GRID_POINTS)
    m_grid = np.linspace(float(log_values.min()) + 1e-4, float(log_values.max()) - 1e-4, M_GRID_POINTS)

    theta = conc[:, None] / (conc[:, None] + 10.0 ** logk[None, :])  # (C, K)
    m = np.log10(baseline_fluorescence + A[None, None, :] * theta[:, :, None])  # (C, K, A)
    # A model mean outside the bins' own fluorescence range contributes nothing to the likelihood,
    # as in the authors' code, where it falls outside the interpolation grid. The range is taken per
    # concentration: with measured gates a sample's top gate can sit well below the global maximum,
    # and testing against the global range instead would make a saturated tight binder infeasible at
    # exactly the low concentrations where the model says it is saturated, inverting the ranking.
    lo_c = (log_values.min(1) + 1e-4)[:, None, None]
    hi_c = (log_values.max(1) - 1e-4)[:, None, None]
    outside = (m < lo_c) | (m > hi_c)
    frac = np.interp(m, m_grid, np.arange(m_grid.size))
    i0 = np.clip(np.floor(frac).astype(int), 0, m_grid.size - 2)
    w = frac - i0
    cidx = np.arange(conc.size)[:, None, None]

    k = library_fraction(reads, a)
    V = reads.shape[0]
    best_k = np.empty(V)
    best_a = np.empty(V)
    sd = np.empty(V)
    for s0 in range(0, V, chunk):
        LL = profile_loglik(reads[s0 : s0 + chunk], a[s0 : s0 + chunk], k[s0 : s0 + chunk], log_values, m_grid)
        for j, Lv in enumerate(LL):
            # An infeasible grid point carries -inf, and 0 * -inf is NaN, which would survive the
            # sum and corrupt the argmax. An interpolation touching an infeasible point is itself
            # infeasible, so it takes -inf.
            with np.errstate(invalid="ignore"):
                val = (1.0 - w) * Lv[cidx, i0] + w * Lv[cidx, i0 + 1]
            val = np.where(np.isnan(val), -np.inf, val)
            total = np.where(outside, 0.0, val).sum(0)  # (K, A)
            ik, ia = np.unravel_index(np.argmax(total), total.shape)
            post = np.exp(total - total[ik, ia])
            pk = post.sum(1) / post.sum()
            best_k[s0 + j] = logk[ik]
            best_a[s0 + j] = A[ia]
            # Spread about the best fit, not the posterior mean: the authors' KD_sigma.
            sd[s0 + j] = float(np.sqrt((pk * (logk - logk[ik]) ** 2).sum()))
    return GridFit(best_k, best_a, sd, logk)


def predicted_mean_bin(
    log10_kd: float,
    amplitude: float,
    concentrations: np.ndarray,
    baseline_fluorescence: float,
    bin_width: float,
    n_bins: int,
) -> np.ndarray:
    """Mean bin (1-based) the fitted curve predicts on the fixed scale, capped at the top bin."""
    theta = concentrations / (concentrations + 10.0**log10_kd)
    m = np.log10(baseline_fluorescence + amplitude * theta) / bin_width
    return 1.0 + np.clip(m, 0.0, n_bins - 1)


def predicted_mean_bin_measured(
    log10_kd: float,
    amplitude: float,
    concentrations: np.ndarray,
    baseline_fluorescence: float,
    log_values: np.ndarray,
    bin_labels: np.ndarray,
) -> np.ndarray:
    """Mean bin the fitted curve predicts when each bin's fluorescence was measured per sample.

    log_values: (C, B) for the same concentrations, ascending across bins. The predicted mean log
    fluorescence is inverted against that concentration's own gate values, so the result stays on
    the observed bin numbering and is comparable with the measured mean bin.
    """
    theta = np.asarray(concentrations, dtype=float) / (np.asarray(concentrations, dtype=float) + 10.0**log10_kd)
    m = np.log10(baseline_fluorescence + amplitude * theta)
    labels = np.asarray(bin_labels, dtype=float)
    return np.array([float(np.interp(m[i], log_values[i], labels)) for i in range(m.size)])
