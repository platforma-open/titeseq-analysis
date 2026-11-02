"""Equation dispatch: the pipeline under "log_fluorescence", and the v2.4.0 path left alone."""

from __future__ import annotations

import json
import math

import numpy as np
import polars as pl
import pytest

from constants import FitParams
from main import main
from pipeline import run

CONCS = [1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6]
BINS = [1, 2, 3, 4]
GATE_WIDTH = 0.65
LOG_LO = math.log10(min(CONCS))
LOG_HI = math.log10(max(CONCS))

LF = FitParams(equation="log_fluorescence")
HILL = FitParams(equation="hill_gate_index")


def _lf_reads(
    clonotype: str,
    true_kd: float,
    *,
    concs: list[float] | None = None,
    baseline: float = 2.0,
    top: float | None = None,
    per_conc: int = 800,
    sd_bins: float = 0.35,
) -> list[dict]:
    """Reads drawn from the model the equation fits: mean log10 fluorescence log10(A*c/(c+K) + B).

    Each concentration's reads are spread over the bins as a Gaussian around the bin index that
    carries the model's fluorescence. Written out here so the test does not reuse src.
    """
    concs = CONCS if concs is None else concs
    span = (len(BINS) - 1) * GATE_WIDTH
    saturated = 10.0 ** (0.85 * span) if top is None else top
    c = np.asarray(concs, dtype=float)
    fluorescence = np.log10((saturated - baseline) * c / (c + true_kd) + baseline)
    targets = fluorescence / GATE_WIDTH + BINS[0]

    rows = []
    for i, conc in enumerate(concs):
        weights = np.exp(-0.5 * ((np.array(BINS, dtype=float) - targets[i]) / sd_bins) ** 2)
        weights /= weights.sum()
        counts = np.round(weights * per_conc).astype(int)
        if counts.sum() == 0:
            counts[0] = 1
        for j, b in enumerate(BINS):
            rows.append(
                {
                    "clonotypeKey": clonotype,
                    "sampleId": f"s_c{i}_b{b}",
                    "concentrationStr": str(conc),
                    "concentration": float(conc),
                    "bin": int(b),
                    "reads": int(counts[j]),
                }
            )
    return rows


def _frame(*row_sets: list[dict]) -> pl.DataFrame:
    return pl.DataFrame([row for rows in row_sets for row in rows])


def _row(outputs: dict[str, pl.DataFrame], clonotype: str) -> dict:
    return outputs["per_clonotype"].filter(pl.col("clonotypeKey") == clonotype).row(0, named=True)


class TestLogFluorescenceFit:
    def test_recovers_a_known_kd_and_rates_it_good(self):
        out = run(_frame(_lf_reads("G1", 1e-8)), params=LF)
        row = _row(out, "G1")
        assert row["affinityClass"] == "Good"
        assert row["fitFailureReason"] is None
        assert math.log10(row["kd"]) == pytest.approx(-8.0, abs=0.3)
        assert row["kdOutOfRange"] is False

    def test_the_hill_coefficient_is_the_fixed_exponent(self):
        row = _row(run(_frame(_lf_reads("G1", 1e-8)), params=LF), "G1")
        assert row["hillCoefficient"] == 1.0
        assert row["hillPlotPosition"] == 1.0

    def test_r2_is_still_reported(self):
        row = _row(run(_frame(_lf_reads("G1", 1e-8)), params=LF), "G1")
        assert row["r2"] is not None
        assert row["r2"] > 0.9

    def test_reports_the_standard_error_that_rated_it_good(self):
        row = _row(run(_frame(_lf_reads("G1", 1e-8)), params=LF), "G1")
        assert row["kdLog10Se"] is not None
        assert math.isfinite(row["kdLog10Se"])
        assert 0.0 < row["kdLog10Se"] <= 0.5

    def test_the_output_columns_are_the_same_under_both_equations(self):
        reads = _frame(_lf_reads("G1", 1e-8))
        lf = run(reads, params=LF)
        hill = run(reads, params=HILL)
        assert set(lf) == set(hill)
        for name in lf:
            assert lf[name].columns == hill[name].columns, name
            assert lf[name].schema == hill[name].schema, name

    def test_the_fitted_curve_comes_back_on_the_bin_index_scale(self):
        out = run(_frame(_lf_reads("G1", 1e-8)), params=LF)
        fitted = out["fitted_mean_bin"]["fittedMeanBin"].to_numpy()
        assert fitted.size == len(CONCS)
        assert fitted.min() >= BINS[0] - 0.5
        assert fitted.max() <= BINS[-1] + 0.5

    def test_a_tighter_binder_reports_a_smaller_kd(self):
        out = run(_frame(_lf_reads("TIGHT", 1e-9), _lf_reads("LOOSE", 1e-7)), params=LF)
        assert _row(out, "TIGHT")["kd"] < _row(out, "LOOSE")["kd"]


class TestBounds:
    def test_a_curve_saturated_across_the_titration_is_pinned_low_and_rated_partial(self):
        # No 0 M point, so nothing in the titration is dim: the curve is flat and bright.
        row = _row(run(_frame(_lf_reads("SAT", 1e-14)), params=LF), "SAT")
        assert row["affinityClass"] == "Partial"
        assert row["fitFailureReason"] is None
        assert row["kdOutOfRange"] is True
        assert row["kdRangeSide"] == "Below range"
        assert math.log10(row["kd"]) == pytest.approx(LOG_LO, abs=1e-9)

    def test_a_pinned_kd_reports_no_standard_error(self):
        # The SE belongs to the unpinned estimate; the reported K_D is a bound.
        row = _row(run(_frame(_lf_reads("SAT", 1e-14)), params=LF), "SAT")
        assert row["kdOutOfRange"] is True
        assert row["kdLog10Se"] is None

    def test_a_flat_dim_curve_is_pinned_high_and_rated_a_non_binder(self):
        row = _row(run(_frame(_lf_reads("DIM", 1e3)), params=LF), "DIM")
        assert (row["affinityClass"], row["fitFailureReason"]) == ("Non-binder", None)
        assert row["kdOutOfRange"] is True
        assert row["kdRangeSide"] == "Above range"
        assert math.log10(row["kd"]) == pytest.approx(LOG_HI, abs=1e-9)

    def test_the_plateau_floor_decides_which_wall_a_flat_mid_bright_curve_is_pinned_to(self):
        # Observed top at 50 percent of the gate range, rising 0.65 decades — under the rise
        # threshold, so the level test alone decides. Saturated at 0.4, a non-binder at 2/3.
        reads = _frame(_lf_reads("MID", 1e-8, top=10 ** (0.50 * 1.95)))
        lenient = _row(run(reads, params=FitParams(equation="log_fluorescence", plateau_floor=0.4)), "MID")
        strict = _row(run(reads, params=FitParams(equation="log_fluorescence", plateau_floor=2 / 3)), "MID")
        assert math.log10(lenient["kd"]) == pytest.approx(LOG_LO, abs=1e-9)
        assert math.log10(strict["kd"]) == pytest.approx(LOG_HI, abs=1e-9)
        assert (lenient["affinityClass"], strict["affinityClass"]) == ("Partial", "Non-binder")
        assert (lenient["kdRangeSide"], strict["kdRangeSide"]) == ("Below range", "Above range")
        assert lenient["kdOutOfRange"] is True and strict["kdOutOfRange"] is True


class TestGates:
    def test_too_few_points_fails_without_fitting(self):
        # Four concentrations, under the equation's fixed minimum of five.
        row = _row(run(_frame(_lf_reads("SHORT", 1e-8, concs=CONCS[:4])), params=LF), "SHORT")
        assert (row["affinityClass"], row["fitFailureReason"]) == ("Failed", "insufficient_points")
        assert row["kd"] is None
        assert row["hillCoefficient"] is None

    def test_the_user_read_floor_is_not_read_by_this_equation(self):
        # min_reads_per_concentration 10_000 would gate every point under v2.4.0's fit.
        reads = _frame(_lf_reads("G1", 1e-8))
        loose = run(reads, params=LF)
        tight = run(reads, params=FitParams(equation="log_fluorescence", min_reads_per_concentration=10_000))
        assert loose["per_clonotype"].equals(tight["per_clonotype"])

    def test_no_bin_mode_is_refused_with_a_reason(self):
        reads = _frame(_lf_reads("G1", 1e-8)).drop("bin")
        with pytest.raises(ValueError, match="bin column"):
            run(reads, params=LF)

    @pytest.mark.parametrize("floor", [0.0, -0.1, 1.5])
    def test_a_plateau_floor_outside_its_range_is_refused(self, floor):
        with pytest.raises(ValueError, match="plateau_floor"):
            run(_frame(_lf_reads("G1", 1e-8)), params=FitParams(equation="log_fluorescence", plateau_floor=floor))


class TestV240PathUnchanged:
    """The v2.4.0 path is what existing users run. Nothing in the "log_fluorescence" path may reach it."""

    def test_an_absent_equation_runs_the_v240_fit(self):
        reads = _frame(_lf_reads("G1", 1e-8))
        assert run(reads)["per_clonotype"].equals(run(reads, params=HILL)["per_clonotype"])

    def test_the_plateau_floor_does_not_reach_the_v240_fit(self):
        reads = _frame(_lf_reads("G1", 1e-8))
        a = run(reads, params=FitParams(equation="hill_gate_index", plateau_floor=2 / 3))
        b = run(reads, params=FitParams(equation="hill_gate_index", plateau_floor=0.05))
        for name in a:
            assert a[name].equals(b[name]), name

    def test_the_v240_fit_reports_no_standard_error(self):
        out = run(_frame(_lf_reads("G1", 1e-8)), params=HILL)
        assert out["per_clonotype"]["kdLog10Se"].null_count() == out["per_clonotype"].height

    def test_the_v240_fit_still_fits_its_own_exponent(self):
        row = _row(run(_frame(_lf_reads("G1", 1e-8)), params=HILL), "G1")
        assert row["hillCoefficient"] != 1.0

    def test_the_two_equations_do_not_agree_on_log_fluorescence_data(self):
        # v2.4.0 fits a linear Hill curve to the mean gate INDEX, which reads systematically tight.
        reads = _frame(_lf_reads("G1", 1e-8))
        lf = math.log10(_row(run(reads, params=LF), "G1")["kd"])
        hill = math.log10(_row(run(reads, params=HILL), "G1")["kd"])
        assert abs(lf - hill) > 0.1


class TestSortFraction:
    """The Adams eq. A3 correction runs before dispatch, so it reaches both equations."""

    def test_the_correction_changes_the_reported_kd(self):
        rows = _lf_reads("G1", 1e-8)
        # Bin 1 takes most of the cells at low concentrations, bin 4 at high ones. Each
        # concentration's four fractions sum to 1, which the validator requires.
        low = dict(zip(BINS, [0.4, 0.3, 0.2, 0.1]))
        high = dict(zip(BINS, [0.1, 0.2, 0.3, 0.4]))
        for row in rows:
            rank = CONCS.index(row["concentration"]) / (len(CONCS) - 1)
            row["sort_fraction"] = (high if rank > 0.5 else low)[row["bin"]]
        reads = pl.DataFrame(rows)
        plain = _row(run(reads.drop("sort_fraction"), params=LF), "G1")
        corrected = _row(run(reads, params=LF, sort_fraction_col="sort_fraction"), "G1")
        assert corrected["kd"] != plain["kd"]
        assert corrected["affinityClass"] in {"Good", "Partial"}


def _with_positions(reads: pl.DataFrame, offset: float = 2.0) -> pl.DataFrame:
    """Measured gate positions at the assumed spacing, shifted to a realistic log10 fluorescence."""
    return reads.with_columns(((pl.col("bin") - 1) * GATE_WIDTH + offset).cast(pl.Float64).alias("bin_fluorescence"))


class TestOutputScale:
    """meanBin and fittedMeanBin are plotted on one axis, so they must share a scale."""

    @pytest.mark.parametrize(
        "measured,lo,hi,gate",
        [
            # Assumed positions: both are the mean gate index, one gate apart per unit.
            (False, BINS[0], BINS[-1], 1.0),
            # Measured positions: both are mean log10 fluorescence, here offset by 2 decades.
            (True, 2.0, 2.0 + (len(BINS) - 1) * GATE_WIDTH, GATE_WIDTH),
        ],
    )
    def test_the_fitted_curve_shares_the_observed_scale(self, measured, lo, hi, gate):
        reads = _frame(_lf_reads("G1", 1e-8))
        out = run(_with_positions(reads) if measured else reads, params=LF)
        joined = out["mean_bin"].join(out["fitted_mean_bin"], on=["clonotypeKey", "concentrationStr"])
        assert joined.height == len(CONCS)
        for name in ("meanBin", "fittedMeanBin"):
            assert joined[name].min() >= lo - 1e-9, name
            assert joined[name].max() <= hi + 1e-9, name
        # Within half a gate at every concentration. Curves on different scales sit gates apart.
        assert (joined["meanBin"] - joined["fittedMeanBin"]).abs().max() < 0.5 * gate


class TestMeasuredGatePositionValidation:
    @pytest.mark.parametrize("bad", [None, math.nan, math.inf])
    def test_a_gate_without_a_finite_position_is_refused(self, bad):
        reads = _with_positions(_frame(_lf_reads("G1", 1e-8)))
        reads = reads.with_columns(
            pl.when(pl.int_range(pl.len()) == 7)
            .then(pl.lit(bad, dtype=pl.Float64))
            .otherwise(pl.col("bin_fluorescence"))
            .alias("bin_fluorescence")
        )
        with pytest.raises(ValueError, match="bin_fluorescence"):
            run(reads, params=LF)

    def test_the_other_equation_does_not_read_the_column(self):
        reads = _with_positions(_frame(_lf_reads("G1", 1e-8)))
        reads = reads.with_columns(pl.lit(None, dtype=pl.Float64).alias("bin_fluorescence"))
        assert run(reads, params=HILL)["per_clonotype"].height == 1


class TestSharedBaseline:
    """With 0 M samples every clonotype is fitted against one baseline: the read-weighted 0 M mean."""

    @staticmethod
    def _zero_molar(clonotype: str, gate: int, reads: int) -> list[dict]:
        # One clonotype alone in its gate at 0 M, so its 0 M signal is exactly that gate's position.
        return [
            {
                "clonotypeKey": clonotype,
                "sampleId": f"s_0M_b{gate}",
                "concentrationStr": "0",
                "concentration": 0.0,
                "bin": gate,
                "reads": reads,
            }
        ]

    def test_every_fit_receives_the_read_weighted_zero_molar_signal(self, monkeypatch):
        import pipeline

        seen: list[float | None] = []
        fit_one = pipeline.fit_one

        def spy(c, y, w, **kwargs):
            seen.append(kwargs["baseline"])
            return fit_one(c, y, w, **kwargs)

        monkeypatch.setattr(pipeline, "fit_one", spy)
        reads = _frame(
            _lf_reads("A", 1e-8),
            _lf_reads("B", 1e-9),
            self._zero_molar("A", gate=1, reads=300),
            self._zero_molar("B", gate=2, reads=100),
        )
        run(reads, params=LF)

        # Gate 1 sits at 0 decades and gate 2 at one gate width. Weighted 300:100 the mean is a
        # quarter of a gate width; unweighted it would be half.
        assert seen == [pytest.approx(0.25 * GATE_WIDTH)] * 2

    def test_without_zero_molar_samples_each_clonotype_fits_its_own(self, monkeypatch):
        import pipeline

        seen: list[float | None] = []
        fit_one = pipeline.fit_one

        def spy(c, y, w, **kwargs):
            seen.append(kwargs["baseline"])
            return fit_one(c, y, w, **kwargs)

        monkeypatch.setattr(pipeline, "fit_one", spy)
        run(_frame(_lf_reads("A", 1e-8), _lf_reads("B", 1e-9)), params=LF)
        assert seen == [None, None]


class TestCliParams:
    """`equation` and `plateau_floor` reach the pipeline through the --params JSON file."""

    def _run(self, tmp_path, params: dict | None) -> pl.DataFrame:
        reads_path = tmp_path / "reads.parquet"
        _frame(_lf_reads("G1", 1e-8)).write_parquet(reads_path)
        pc_path = tmp_path / "pc.parquet"
        argv = [
            "--reads",
            str(reads_path),
            "--out-per-clonotype",
            str(pc_path),
            "--out-mean-bin",
            str(tmp_path / "mb.parquet"),
            "--out-fitted-mean-bin",
            str(tmp_path / "fmb.parquet"),
            "--out-concentration-value",
            str(tmp_path / "cv.parquet"),
        ]
        if params is not None:
            params_path = tmp_path / "params.json"
            params_path.write_text(json.dumps(params))
            argv += ["--params", str(params_path)]
        assert main(argv) == 0
        return pl.read_parquet(pc_path)

    def test_log_fluorescence_arrives_through_the_params_file(self, tmp_path):
        frame = self._run(tmp_path, {"equation": "log_fluorescence", "plateau_floor": 2 / 3})
        assert frame.row(0, named=True)["hillCoefficient"] == 1.0

    def test_an_absent_equation_runs_the_v240_fit(self, tmp_path):
        # A pre-existing project's args carry no `equation`. It must not move to "log_fluorescence".
        frame = self._run(tmp_path, {"min_reads_per_concentration": 3, "min_concentration_points": 5})
        assert frame.row(0, named=True)["hillCoefficient"] != 1.0

    def test_no_params_file_at_all_runs_the_v240_fit(self, tmp_path):
        assert self._run(tmp_path, None).row(0, named=True)["hillCoefficient"] != 1.0

    def test_an_absent_plateau_floor_is_harmless(self, tmp_path):
        frame = self._run(tmp_path, {"equation": "hill_gate_index", "n_min": 0.5, "n_max": 2.0})
        assert frame.height == 1

    def test_an_unknown_equation_is_refused_by_name(self, tmp_path):
        with pytest.raises(ValueError, match="unknown equation"):
            self._run(tmp_path, {"equation": "poisson"})


class TestLinearGateFluorescence:
    """A linear MFI where log10 is required moves K_D by decades, or overflows the fit."""

    @pytest.mark.parametrize("offset", [0.0, 2.0])
    def test_linear_fluorescence_is_refused_and_named(self, offset):
        # 10^((bin - 1) * 0.65 + offset): 1 to 89 a.u., and 100 to 8,900 a.u.
        reads = _with_positions(_frame(_lf_reads("G1", 1e-8)), offset=offset)
        reads = reads.with_columns((10.0 ** pl.col("bin_fluorescence")).alias("bin_fluorescence"))
        with pytest.raises(ValueError, match="Gate fluorescence column .* linear fluorescence"):
            run(reads, params=LF)

    def test_log10_fluorescence_up_to_the_limit_is_accepted(self):
        # Gates at 5.05 to 7.0 log10, the top one on the limit.
        reads = _with_positions(_frame(_lf_reads("G1", 1e-8)), offset=7.0 - 3 * GATE_WIDTH)
        assert _row(run(reads, params=LF), "G1")["affinityClass"] in {"Good", "Partial"}


class TestR2Gate:
    def test_off_by_default(self):
        assert _row(run(_frame(_lf_reads("G1", 1e-8)), params=LF), "G1")["affinityClass"] == "Good"

    def test_an_in_range_fit_below_the_gate_fails_low_r2(self):
        # A sampled curve never fits perfectly, so a gate at 1.0 fails every in-range fit.
        gated = FitParams(equation="log_fluorescence", r2_gate_min=1.0)
        row = _row(run(_frame(_lf_reads("G1", 1e-8)), params=gated), "G1")
        assert (row["affinityClass"], row["fitFailureReason"]) == ("Failed", "low_r2")
        assert row["kdRangeSide"] == "In range"
        assert row["r2"] is not None and row["r2"] < 1.0

    def test_a_bound_is_not_gated(self):
        gated = FitParams(equation="log_fluorescence", r2_gate_min=1.0)
        row = _row(run(_frame(_lf_reads("SAT", 1e-14)), params=gated), "SAT")
        assert row["affinityClass"] == "Partial"

    def test_the_mean_bin_fit_ignores_the_gate(self):
        gated = FitParams(equation="hill_gate_index", r2_gate_min=1.0)
        assert _row(run(_frame(_lf_reads("G1", 1e-8)), params=gated), "G1")["affinityClass"] == "Good"


class TestRangeSide:
    def test_an_in_range_fit_reports_in_range(self):
        assert _row(run(_frame(_lf_reads("G1", 1e-8)), params=LF), "G1")["kdRangeSide"] == "In range"

    def test_the_mean_bin_fit_writes_no_side(self):
        # "hill_gate_index" keeps the v2.4.0 columns; kdOutOfRange carries its range fact.
        assert _row(run(_frame(_lf_reads("G1", 1e-8)), params=HILL), "G1")["kdRangeSide"] is None

    def test_a_row_without_a_kd_has_no_side(self):
        row = _row(run(_frame(_lf_reads("SHORT", 1e-8, concs=CONCS[:4])), params=LF), "SHORT")
        assert row["kd"] is None and row["kdRangeSide"] is None


class TestGridWarningReadsTheFilteredFrame:
    def test_another_antigens_samples_do_not_count(self, capsys):
        # G1 sits in every sample of antigen X. Antigen Y's samples hold only OTHER, so before the
        # filter G1 would miss Y's whole grid.
        x = [dict(r, antigen="X") for r in _lf_reads("G1", 1e-8)]
        y = [
            dict(r, antigen="Y", sampleId="y_" + r["sampleId"], concentrationStr="5e-9", concentration=5e-9)
            for r in _lf_reads("OTHER", 1e-8)[:4]
        ]
        run(_frame(x, y), params=LF, target_antigen="X", antigen_column_ref="antigen")
        assert "(bin, concentration) samples" not in capsys.readouterr().out
