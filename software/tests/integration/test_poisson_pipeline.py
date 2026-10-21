"""Behavioral tests for the default bin-mode model (poisson_log) through pipeline.run, on read-level
simulations built like Adams et al. 2016 Appendix G (see tests/conftest.py)."""

from __future__ import annotations

import io
from contextlib import redirect_stdout

import numpy as np
import polars as pl
import pytest

from constants import FitParams
from io_layer import InputValidationError
from pipeline import run


def _run(reads, **kw):
    with redirect_stdout(io.StringIO()):
        return run(reads, **kw)


@pytest.fixture(scope="module")
def fitted(simulated_library):
    out = _run(simulated_library.reads, params=FitParams(bin_width_decades=0.7), cell_count_col="cells")
    pc = out["per_clonotype"].with_columns(pl.col("kd").log10().alias("log10_kd"))
    return out, {row["clonotypeKey"]: row for row in pc.iter_rows(named=True)}


def test_recovers_kd_across_the_titration_range(simulated_library, fitted):
    # Eq. E2 equates the bin-averaged value with mean log fluorescence; cells spilling past the
    # outer bins make weak binders read up to ~0.2 decades weaker. That is the method, not noise.
    _, rows = fitted
    for key, true in simulated_library.true_log10_kd.items():
        if key.startswith("B"):
            assert rows[key]["log10_kd"] == pytest.approx(true, abs=0.25), key


def test_kd_does_not_depend_on_expression(fitted):
    # The point of titrating: two variants with the same K_D and 2x different display give the same K_D.
    _, rows = fitted
    for j in range(7):
        assert abs(rows[f"B{j}0"]["log10_kd"] - rows[f"B{j}1"]["log10_kd"]) < 0.2


def test_in_range_binders_are_good_with_tight_uncertainty(fitted):
    _, rows = fitted
    for key, row in rows.items():
        if key.startswith("B"):
            assert row["affinityClass"] == "Good", key
            assert row["kdCensor"] is None
            assert row["kdLog10Sd"] < 0.3


def test_binder_tighter_than_the_lowest_concentration_is_a_lower_bound(fitted):
    _, rows = fitted
    tight = rows["TIGHT"]
    assert tight["kdCensor"] == "below"
    assert tight["log10_kd"] == pytest.approx(-10.0)
    assert tight["affinityClass"] == "Good"


def test_flat_curve_is_reported_as_above_the_range_not_as_a_binder(fitted):
    # Without the amplitude prior a flat curve fits as a very tight binder (Adams 2016 deviation D5).
    _, rows = fitted
    nb = rows["NONBINDER"]
    assert nb["kdCensor"] == "above"
    assert nb["log10_kd"] == pytest.approx(-5.5)
    assert (nb["affinityClass"], nb["fitFailureReason"]) == ("Failed", "kd_above_range")


def test_fitted_curve_is_emitted_only_for_non_failed_and_stays_within_bins(fitted):
    out, rows = fitted
    fmb = out["fitted_mean_bin"]
    assert "NONBINDER" not in set(fmb["clonotypeKey"].to_list())
    assert {"B00", "TIGHT"} <= set(fmb["clonotypeKey"].to_list())
    assert fmb["fittedMeanBin"].min() >= 1.0 and fmb["fittedMeanBin"].max() <= 4.0


def test_sort_fractions_give_the_same_kd_when_reads_do_not_exceed_cells(simulate):
    sim = simulate({"A": (-8.0, 60.0), "B": (-7.0, 60.0), "N": (-2.0, 60.0)}, reads_per_sample=500, seed=3)
    # Fractions are per sample: compute them from one row per sample, not per clonotype row.
    per_sample = sim.reads.select("sampleId", "concentrationStr", "cells").unique()
    per_sample = per_sample.with_columns(
        (pl.col("cells") / pl.col("cells").sum().over("concentrationStr")).alias("frac")
    )
    reads = sim.reads.join(per_sample.select("sampleId", "frac"), on="sampleId")
    params = FitParams(bin_width_decades=0.7, min_concentration_points=3)
    by_cells = _run(reads, params=params, cell_count_col="cells")["per_clonotype"]
    by_frac = _run(reads.drop("cells"), params=params, sort_fraction_col="frac")["per_clonotype"]
    assert by_cells["kd"].to_list() == pytest.approx(by_frac["kd"].to_list())


class TestInputsTheModelNeeds:
    @pytest.fixture(scope="class")
    def small(self, simulate):
        return simulate({"A": (-8.0, 60.0), "N": (-2.0, 60.0)}, seed=5).reads

    def test_bin_mode_without_cell_counts_is_rejected(self, small):
        with pytest.raises(InputValidationError, match="cells sorted into each sample"):
            _run(small.drop("cells"))

    def test_uncorrected_runs_only_when_asked_and_says_so(self, small, capsys):
        out = run(small.drop("cells"), params=FitParams(allow_uncorrected=True))
        assert out["per_clonotype"].height == 2
        assert "uncorrected" in capsys.readouterr().out

    def test_missing_zero_control_is_rejected(self, small):
        with pytest.raises(InputValidationError, match="0 M"):
            _run(small.filter(pl.col("concentration") > 0), cell_count_col="cells")

    def test_baseline_override_replaces_the_zero_control(self, small):
        params = FitParams(bin_width_decades=0.7, baseline_fluorescence=2.5)
        out = _run(small.filter(pl.col("concentration") > 0), params=params, cell_count_col="cells")
        a = out["per_clonotype"].filter(pl.col("clonotypeKey") == "A")
        assert np.log10(a["kd"][0]) == pytest.approx(-8.0, abs=0.25)

    def test_plateau_prior_below_baseline_is_rejected(self, small):
        with pytest.raises(InputValidationError, match="min_plateau_fraction"):
            _run(small, params=FitParams(min_plateau_fraction=0.05), cell_count_col="cells")

    @pytest.mark.parametrize(
        "corrupt, message",
        [
            (lambda df: df.with_columns(pl.lit(-1.0).alias("cells")), "finite and >= 0"),
            (
                lambda df: df.with_columns(
                    pl.when(pl.col("clonotypeKey") == "A").then(1.0).otherwise(pl.col("cells")).alias("cells")
                ),
                "one value per sample",
            ),
            (lambda df: df.with_columns(pl.lit(None, dtype=pl.Float64).alias("cells")), "null"),
        ],
        ids=["negative", "differs-within-sample", "null"],
    )
    def test_invalid_cell_counts_are_rejected(self, small, corrupt, message):
        with pytest.raises(InputValidationError, match=message):
            _run(corrupt(small), cell_count_col="cells")


@pytest.mark.parametrize(
    "params, message",
    [
        (FitParams(bin_width_decades=0.0), "bin_width_decades"),
        (FitParams(min_plateau_fraction=1.0), "min_plateau_fraction"),
    ],
    ids=["zero-bin-width", "plateau-fraction-not-below-1"],
)
def test_invalid_model_settings_are_rejected(simulate, params, message):
    reads = simulate({"A": (-8.0, 60.0)}, seed=7).reads
    with pytest.raises(InputValidationError, match=message):
        _run(reads, params=params, cell_count_col="cells")


def test_uncertain_in_range_kd_is_partial(simulated_library):
    # kd_sd_good is the Good/Partial line: with a zero tolerance every in-range K_D is Partial,
    # while bounds keep their own classes.
    out = _run(simulated_library.reads, params=FitParams(bin_width_decades=0.7, kd_sd_good=0.0), cell_count_col="cells")
    rows = {r["clonotypeKey"]: r for r in out["per_clonotype"].iter_rows(named=True)}
    assert rows["B30"]["affinityClass"] == "Partial"
    assert rows["TIGHT"]["affinityClass"] == "Good"
    assert rows["NONBINDER"]["affinityClass"] == "Failed"


def test_clonotype_seen_at_too_few_concentrations_is_not_fitted(simulate):
    sim = simulate({"A": (-8.0, 60.0), "N": (-2.0, 60.0)}, seed=9)
    concs = sorted(sim.reads["concentration"].unique().to_list())
    # Keep "A" at 0 M and 3 nonzero concentrations only: below min_concentration_points (5).
    keep = pl.col("clonotypeKey").ne("A") | pl.col("concentration").is_in(concs[:4])
    out = _run(sim.reads.filter(keep), params=FitParams(bin_width_decades=0.7), cell_count_col="cells")
    a = out["per_clonotype"].filter(pl.col("clonotypeKey") == "A")
    assert (a["affinityClass"][0], a["fitFailureReason"][0]) == ("Failed", "insufficient_points")
    assert a["kd"][0] is None and a["kdCensor"][0] is None
    assert "A" not in set(out["fitted_mean_bin"]["clonotypeKey"].to_list())
