"""Golden tests: the block reproduces Adams et al. 2016 (eLife 5:e23156) from the authors' read counts.

Each replicate runs through `pipeline.run` with the experiment's bin width and the authors' amplitude
prior, as an operator would set them. Assertions compare to the paper (Table I, WT K_D, Fig 4C) and to
the authors' published per-variant K_D. Downloads the data on first run (see conftest.py).
"""

from __future__ import annotations

import io
from contextlib import redirect_stdout

import numpy as np
import polars as pl
import pytest
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
from pipeline import run

pytestmark = [pytest.mark.golden, pytest.mark.slow]


def _fit(rep, **settings) -> pl.DataFrame:
    params = FitParams(bin_width_decades=rep.bin_width_decades, **settings)
    with redirect_stdout(io.StringIO()):
        out = run(rep.reads, params=params, cell_count_col="cell_count")
    pc = out["per_clonotype"].with_columns(pl.col("kd").log10().alias("log10_kd"))
    return rep.published.join(pc, on="clonotypeKey", how="left")


@pytest.fixture(scope="module")
def paper_settings(adams_replicates):
    """The reproduction proper: the authors' amplitude prior, and every variant fitted (no read gate),
    as Appendix E does."""
    return {
        rep.number: _fit(
            rep,
            min_plateau_fraction=rep.min_plateau_fraction,
            min_reads_per_concentration=1,
            min_concentration_points=1,
        )
        for rep in adams_replicates
    }


@pytest.fixture(scope="module")
def default_settings(adams_replicates):
    """Block defaults throughout; only the bin width comes from the experiment."""
    return {rep.number: _fit(rep) for rep in adams_replicates}


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


@pytest.mark.parametrize("number", [1, 2, 3])
def test_per_variant_kd_matches_published(paper_settings, number):
    # Replicates of the same experiment agree with each other at r ~ 0.88; matching the authors' table
    # well above that means the same method, not a similar one. The residual is the authors'
    # Nelder-Mead profile and grid interpolation, which this implementation solves exactly.
    j = paper_settings[number].filter(pl.col("published_log10_kd").is_not_nan() & pl.col("log10_kd").is_not_null())
    ours, pub = j["log10_kd"].to_numpy(), j["published_log10_kd"].to_numpy()
    diff = np.abs(ours - pub)
    assert j.height >= 0.99 * paper_settings[number].height
    assert np.corrcoef(ours, pub)[0, 1] >= 0.95
    assert np.mean(diff < 0.25) >= 0.92
    assert np.median(diff) <= 0.02


def test_wt_kd_matches_paper(paper_settings):
    # Results D: WT averaged over its synonymous variants and the three replicates is 10^-8.87.
    wt = pl.concat(
        f.filter((pl.col("cdr1h_aa") == WT_CDR1H_AA) & (pl.col("cdr3h_aa") == WT_CDR3H_AA))
        for f in paper_settings.values()
    )
    assert wt.height >= 90
    assert wt["log10_kd"].mean() == pytest.approx(WT_LOG10_KD, abs=0.1)


def test_table_i_clones_match_paper(paper_settings, adams_clone_keys):
    # Table I Tite-Seq column: each clone's replicate mean, within 0.3 decades or two standard errors.
    means = _clone_means(paper_settings, adams_clone_keys)
    misses = [
        (c.name, round(means[c.name], 2), c.tite_seq_log10_kd)
        for c in TABLE_I
        if not abs(means[c.name] - c.tite_seq_log10_kd) <= max(0.3, 2 * (c.tite_seq_se or 0.0))
    ]
    assert misses == []


def test_table_i_agrees_with_flow_like_fig_4c(paper_settings, adams_clone_keys):
    # Fig 4C: R^2 = 0.89 between Tite-Seq and flow K_D over the Table I clones, bounds at the limits.
    # The 0.04 allowance is clone C93, replicate 2: the authors report 10^-5.9 where their other two
    # replicates and this fit give 10^-8.3 to 10^-8.5, the source of Table I's +-0.70 for C93.
    assert _r2(_clone_means(paper_settings, adams_clone_keys)) >= FIG4C_R2 - 0.04


def test_non_binders_are_reported_as_upper_bounds(paper_settings):
    # Variants the authors put at >= 10^-5 M are bounds, not measurements; the block must say so.
    for fits in paper_settings.values():
        nb = fits.filter(pl.col("published_log10_kd") >= LOG10_KD_UPPER_LIMIT)
        assert (nb["kdCensor"] == "above").mean() >= 0.95


def test_reported_kd_stays_within_titration_range(paper_settings):
    for fits in paper_settings.values():
        kd = fits["log10_kd"].drop_nulls()
        assert kd.min() >= LOG10_KD_LOWER_LIMIT - 1e-9
        assert kd.max() <= LOG10_KD_UPPER_LIMIT + 1e-9


def test_default_settings_still_track_flow(default_settings, adams_clone_keys):
    # Block defaults (q = 0.65 instead of the authors' 0.63-0.70; read gate on): Table I must still
    # agree with flow, and every clone must get a K_D or a bound.
    means = _clone_means(default_settings, adams_clone_keys)
    assert all(np.isfinite(v) for v in means.values())
    assert _r2(means) >= 0.80
