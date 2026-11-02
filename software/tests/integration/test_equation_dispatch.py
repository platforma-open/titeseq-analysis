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
        assert math.log10(row["kd"]) == pytest.approx(LOG_LO, abs=1e-9)

    def test_a_flat_dim_curve_is_pinned_high_and_rated_partial(self):
        row = _row(run(_frame(_lf_reads("DIM", 1e3)), params=LF), "DIM")
        assert row["affinityClass"] == "Partial"
        assert row["kdOutOfRange"] is True
        assert math.log10(row["kd"]) == pytest.approx(LOG_HI, abs=1e-9)

    def test_the_plateau_floor_decides_which_wall_a_flat_mid_bright_curve_is_pinned_to(self):
        # Observed top at 50 percent of the gate range, rising 0.65 decades — under the rise
        # threshold, so the level test alone decides. Saturated at 0.4, a non-binder at 2/3.
        reads = _frame(_lf_reads("MID", 1e-8, top=10 ** (0.50 * 1.95)))
        lenient = _row(run(reads, params=FitParams(equation="log_fluorescence", plateau_floor=0.4)), "MID")
        strict = _row(run(reads, params=FitParams(equation="log_fluorescence", plateau_floor=2 / 3)), "MID")
        assert math.log10(lenient["kd"]) == pytest.approx(LOG_LO, abs=1e-9)
        assert math.log10(strict["kd"]) == pytest.approx(LOG_HI, abs=1e-9)
        assert lenient["affinityClass"] == strict["affinityClass"] == "Partial"
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
    """The v2.4.0 path is what existing users run. Nothing in this change may reach it."""

    def test_an_absent_equation_runs_the_v240_fit(self):
        reads = _frame(_lf_reads("G1", 1e-8))
        assert run(reads)["per_clonotype"].equals(run(reads, params=HILL)["per_clonotype"])

    def test_the_plateau_floor_does_not_reach_the_v240_fit(self):
        reads = _frame(_lf_reads("G1", 1e-8))
        a = run(reads, params=FitParams(equation="hill_gate_index", plateau_floor=2 / 3))
        b = run(reads, params=FitParams(equation="hill_gate_index", plateau_floor=0.05))
        for name in a:
            assert a[name].equals(b[name]), name

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

    def test_the_new_equation_arrives_through_the_params_file(self, tmp_path):
        frame = self._run(tmp_path, {"equation": "log_fluorescence", "plateau_floor": 2 / 3})
        assert frame.row(0, named=True)["hillCoefficient"] == 1.0

    def test_an_absent_equation_runs_the_v240_fit(self, tmp_path):
        # A pre-existing project's args carry no `equation`. It must not move to the new fit.
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
