"""Golden tests: the block reproduces Adams et al. 2016 (eLife 5:e23156) from the authors' read counts.

Each replicate runs through `pipeline.run` with the experiment's bin width and the authors' amplitude
prior, as an operator would set them. Assertions compare to the paper (Table I, WT K_D, Fig 4C) and to
the authors' published per-variant K_D. Downloads the data on first run (see conftest.py).

Two tiers run the same assertions:
  fast  ~250 variants per replicate (Table I clones, WT synonyms, a stratified sample); the rest of
        the library is pooled so read totals stay exact, and B is the authors' value.
        `uv run pytest -m "golden and not slow"`
  full  every variant, B computed by the block. `uv run pytest -m golden`
Fits are cached on disk, keyed by the source code and settings, so unchanged code reruns in seconds.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from golden_data import GoldenReplicate, fast_subset, fit_golden
from paper_adams2016 import (
    FIG4C_R2,
    LOG10_KD_LOWER_LIMIT,
    LOG10_KD_UPPER_LIMIT,
    TABLE_I,
    WT_CDR1H_AA,
    WT_CDR3H_AA,
    WT_LOG10_KD,
)

from constants import FitParams
from poisson_fit import baseline, exposure
from poisson_pipeline import build_tensors

pytestmark = pytest.mark.golden

TIERS = [pytest.param("fast", id="fast"), pytest.param("full", id="full", marks=pytest.mark.slow)]


def _paper_params(rep: GoldenReplicate, *, baseline_fluorescence: float | None = None) -> FitParams:
    """The authors' settings: their amplitude prior, and every variant fitted (no read gate)."""
    return FitParams(
        allow_computed_baseline=True,
        bin_width_decades=rep.bin_width_decades,
        min_plateau_fraction=rep.min_plateau_fraction,
        min_reads_per_concentration=1,
        min_concentration_points=1,
        baseline_fluorescence=baseline_fluorescence,
    )


@pytest.fixture(scope="module")
def fits_full(adams_replicates: list[GoldenReplicate]) -> dict[int, pl.DataFrame]:
    return {rep.number: fit_golden(rep, _paper_params(rep)) for rep in adams_replicates}


@pytest.fixture(scope="module")
def fits_fast(adams_replicates: list[GoldenReplicate], adams_clone_keys: dict[str, str]) -> dict[int, pl.DataFrame]:
    return {
        rep.number: fit_golden(
            rep,
            _paper_params(rep, baseline_fluorescence=rep.baseline_fluorescence),
            fast_subset(rep, adams_clone_keys),
        )
        for rep in adams_replicates
    }


@pytest.fixture(params=TIERS)
def fits(request: pytest.FixtureRequest) -> dict[int, pl.DataFrame]:
    return request.getfixturevalue(f"fits_{request.param}")


def _clone_means(fits: dict[int, pl.DataFrame], keys: dict[str, str]) -> dict[str, float]:
    means = {}
    for clone in TABLE_I:
        vals = [f.filter(pl.col("clonotypeKey") == keys[clone.name])["log10_kd"] for f in fits.values()]
        vals = [v[0] for v in vals if v.len() and v[0] is not None]
        means[clone.name] = float(np.mean(vals)) if vals else float("nan")
    return means


def _r2(pred: dict[str, float]) -> float:
    """Squared Pearson correlation with flow, bounds at the limits: Table I's own values give 0.90,
    the 0.89 printed in Fig 4C. Clones without a value count as missing, so they fail the assertion."""
    flow = np.array([c.flow_log10_kd for c in TABLE_I])
    est = np.array([pred[c.name] for c in TABLE_I])
    return float(np.corrcoef(est, flow)[0, 1] ** 2)


# (r, share within 0.25 decades, median abs difference). The fast sample holds equal counts per
# published-K_D stratum, so it over-weights the mid-range variants where the authors' own replicates
# disagree most; the same fits score r 0.93-0.99 there and 0.95-0.99 over the whole library.
AGREEMENT = {"full": (0.95, 0.92, 0.02), "fast": (0.92, 0.80, 0.05)}


@pytest.mark.parametrize("number", [1, 2, 3])
def test_per_variant_kd_matches_published(fits, number, request):
    # Replicates of the same experiment agree with each other at r ~ 0.88; matching the authors' table
    # well above that means the same method, not a similar one. The residual is the authors'
    # Nelder-Mead profile and grid interpolation, which this implementation solves exactly.
    min_r, min_close, max_median = AGREEMENT[request.node.callspec.params["fits"]]
    f = fits[number]
    j = f.filter(pl.col("published_log10_kd").is_not_nan() & pl.col("log10_kd").is_not_null())
    ours, pub = j["log10_kd"].to_numpy(), j["published_log10_kd"].to_numpy()
    diff = np.abs(ours - pub)
    assert j.height >= 0.99 * f.filter(pl.col("published_log10_kd").is_not_nan()).height
    assert np.corrcoef(ours, pub)[0, 1] >= min_r
    assert np.mean(diff < 0.25) >= min_close
    assert np.median(diff) <= max_median


def test_wt_kd_matches_paper(fits):
    # Results D: WT averaged over its synonymous variants and the three replicates is 10^-8.87.
    wt = pl.concat(
        f.filter((pl.col("cdr1h_aa") == WT_CDR1H_AA) & (pl.col("cdr3h_aa") == WT_CDR3H_AA)) for f in fits.values()
    )
    assert wt.height >= 90
    assert wt["log10_kd"].mean() == pytest.approx(WT_LOG10_KD, abs=0.1)


def test_table_i_clones_match_paper(fits, adams_clone_keys):
    # Table I Tite-Seq column: each clone's replicate mean, within 0.3 decades or two standard errors.
    means = _clone_means(fits, adams_clone_keys)
    misses = [
        (c.name, round(means[c.name], 2), c.tite_seq_log10_kd)
        for c in TABLE_I
        if not abs(means[c.name] - c.tite_seq_log10_kd) <= max(0.3, 2 * (c.tite_seq_se or 0.0))
    ]
    assert misses == []


def test_table_i_agrees_with_flow_like_fig_4c(fits, adams_clone_keys):
    # Fig 4C: R^2 = 0.89 between Tite-Seq and flow K_D over the Table I clones, bounds at the limits.
    # The 0.04 allowance is clone C93, replicate 2: the authors report 10^-5.9 where their other two
    # replicates and this fit give 10^-8.3 to 10^-8.5, the source of Table I's +-0.70 for C93.
    assert _r2(_clone_means(fits, adams_clone_keys)) >= FIG4C_R2 - 0.04


def test_non_binders_are_reported_as_upper_bounds(fits):
    # Variants the authors put at >= 10^-5 M are bounds, not measurements; the block must say so.
    for f in fits.values():
        nb = f.filter(pl.col("published_log10_kd") >= LOG10_KD_UPPER_LIMIT)
        assert (nb["kdCensor"] == "above").mean() >= 0.95


def test_reported_kd_stays_within_titration_range(fits):
    for f in fits.values():
        kd = f["log10_kd"].drop_nulls()
        assert kd.min() >= LOG10_KD_LOWER_LIMIT - 1e-9
        assert kd.max() <= LOG10_KD_UPPER_LIMIT + 1e-9


@pytest.mark.parametrize("number", [1, 2, 3])
def test_baseline_from_zero_m_reads_matches_the_authors(adams_replicates, number):
    # Appendix E: B is the mean fluorescence of cells at 0 M. The block's estimate from the 0 M reads
    # must equal the authors' value, which the fast tier passes in directly.
    rep = adams_replicates[number - 1]
    t = build_tensors(rep.reads, "cell_count")
    a, h = exposure(t.reads.sum(0), t.sample_cells)
    log_values = rep.bin_width_decades * (t.bin_labels - 1.0)
    b = baseline((t.reads * h[None])[:, 0, :], a[0], log_values)
    assert b == pytest.approx(rep.baseline_fluorescence, abs=2e-3)


@pytest.mark.slow
def test_default_settings_still_track_flow(adams_replicates, adams_clone_keys):
    # Block defaults (q = 0.65 instead of the authors' 0.63-0.70; read gate on): Table I must still
    # agree with flow, and every clone must get a K_D or a bound.
    fits = {
        rep.number: fit_golden(rep, FitParams(allow_computed_baseline=True, bin_width_decades=rep.bin_width_decades))
        for rep in adams_replicates
    }
    means = _clone_means(fits, adams_clone_keys)
    assert all(np.isfinite(v) for v in means.values())
    assert _r2(means) >= 0.80
