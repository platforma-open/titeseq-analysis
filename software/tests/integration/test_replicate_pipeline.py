"""Sample replicates through the whole pipeline, on a synthetic library with known K_D.

Each replicate is a full titration of its own: its own samples, its own gate depths and its own
sort fractions. The expectations are the mean of the per-replicate fits, the counts, and — for a
run with no replicate column — that nothing moves at all.
"""

from __future__ import annotations

import json
import math

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from constants import FitParams
from io_layer import InputValidationError
from main import main
from pipeline import run
from replicate_combine import Filters, ReplicateInputs

CONCS = [1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6]
BINS = [1, 2, 3, 4]
GATE_WIDTH = 0.65
LOG_LO = math.log10(min(CONCS))
LOG_HI = math.log10(max(CONCS))

LF = FitParams(equation="log_fluorescence")
HILL = FitParams(equation="hill_gate_index")


def _reads(nt: str, true_kd: float, replicate: str, *, per_conc: int = 800, sd_bins: float = 0.35) -> list[dict]:
    """Reads drawn from log10(A*c/(c+K) + B), spread over the gates as a Gaussian in gate index."""
    span = (len(BINS) - 1) * GATE_WIDTH
    saturated, baseline = 10.0 ** (0.85 * span), 2.0
    c = np.asarray(CONCS, dtype=float)
    targets = np.log10((saturated - baseline) * c / (c + true_kd) + baseline) / GATE_WIDTH + BINS[0]
    rows = []
    for i, conc in enumerate(CONCS):
        weights = np.exp(-0.5 * ((np.array(BINS, dtype=float) - targets[i]) / sd_bins) ** 2)
        counts = np.round(weights / weights.sum() * per_conc).astype(int)
        for j, b in enumerate(BINS):
            rows.append(
                {
                    "clonotypeKey": nt,
                    "sampleId": f"{replicate}_c{i}_b{b}",
                    "concentrationStr": str(conc),
                    "concentration": float(conc),
                    "bin": int(b),
                    "reads": int(counts[j]),
                    "replicate": replicate,
                    "sort_fraction": 1.0 / len(BINS),
                }
            )
    return rows


def _frame(design: dict[str, dict[str, float]], **depth: int) -> pl.DataFrame:
    """`design` maps a replicate to the true K_D of each variant in it."""
    rows = [
        row
        for replicate, variants in design.items()
        for nt, kd in variants.items()
        for row in _reads(nt, kd, replicate, per_conc=depth.get(replicate, 800))
    ]
    return pl.DataFrame(rows)


THREE = {
    "A": {"v1": 1e-8, "v2": 3e-9},
    "B": {"v1": 2e-8, "v2": 3e-9},
    "C": {"v1": 4e-8, "v2": 3e-9},
}

# The same library sorted three times: every replicate must produce the same fit.
IDENTICAL = {r: {"v1": 1e-8, "v2": 3e-9} for r in ("A", "B", "C")}

REPS = ReplicateInputs(column="replicate")


def _kd(frame: pl.DataFrame, key: str, column: str = "kd") -> float | None:
    return frame.filter(pl.col("clonotypeKey") == key)[column].item()


class TestNoReplicateColumn:
    """With no replicate column asked for, every output is the one the block always produced."""

    def test_the_output_set_and_schemas_do_not_move(self):
        outputs = run(_frame({"A": {"v1": 1e-8}}), params=LF)
        assert set(outputs) == {"per_clonotype", "mean_bin", "fitted_mean_bin", "concentration_value"}
        assert outputs["mean_bin"].columns == ["clonotypeKey", "concentrationStr", "meanBin"]
        assert outputs["fitted_mean_bin"].columns == ["clonotypeKey", "concentrationStr", "fittedMeanBin"]
        assert "replicate" not in outputs["per_clonotype"].columns
        assert "kdLog10Sd" not in outputs["per_clonotype"].columns
        assert not any(c.startswith("replicate") for c in outputs["per_clonotype"].columns)

    def test_a_replicate_column_in_the_reads_changes_nothing_on_its_own(self):
        """The column reaches the frame whenever the workflow joins it; only the flag may act."""
        reads = _frame(THREE)
        with_column = run(reads, params=LF)
        without_column = run(reads.drop("replicate"), params=LF)
        for name in ("per_clonotype", "mean_bin", "fitted_mean_bin", "concentration_value"):
            assert_frame_equal(with_column[name], without_column[name])

    def test_the_cli_writes_only_the_legacy_frames(self, tmp_path):
        _frame({"A": {"v1": 1e-8}}).write_parquet(tmp_path / "reads.parquet")
        out = tmp_path / "out"
        out.mkdir()
        rc = main(
            [
                "--reads", str(tmp_path / "reads.parquet"),
                "--out-per-clonotype", str(out / "per_clonotype.tsv"),
                "--out-mean-bin", str(out / "mean_bin.tsv"),
                "--out-fitted-mean-bin", str(out / "fitted_mean_bin.tsv"),
                "--out-concentration-value", str(out / "concentration_value.tsv"),
            ]
        )  # fmt: skip
        assert rc == 0
        assert sorted(p.name for p in out.iterdir()) == [
            "concentration_value.tsv",
            "fitted_mean_bin.tsv",
            "mean_bin.tsv",
            "per_clonotype.tsv",
        ]

    def test_the_cli_refuses_a_replicate_setting_without_the_column(self):
        with pytest.raises(SystemExit):
            main(
                [
                    "--reads", "r.parquet",
                    "--out-per-clonotype", "a", "--out-mean-bin", "b",
                    "--out-fitted-mean-bin", "c", "--out-concentration-value", "d",
                    "--replicate-error-model", "phillips",
                ]
            )  # fmt: skip

    def test_the_cli_refuses_the_column_without_its_output(self):
        with pytest.raises(SystemExit):
            main(
                [
                    "--reads", "r.parquet",
                    "--out-per-clonotype", "a", "--out-mean-bin", "b",
                    "--out-fitted-mean-bin", "c", "--out-concentration-value", "d",
                    "--replicate-column", "replicate",
                ]
            )  # fmt: skip


class TestPerReplicateFitting:
    def test_each_replicate_is_fitted_on_its_own(self):
        outputs = run(_frame(THREE), params=LF, replicates=REPS)
        per_replicate = outputs["per_replicate"]
        assert per_replicate.height == 6
        v1 = per_replicate.filter(pl.col("clonotypeKey") == "v1").sort("replicate")
        assert v1["replicate"].to_list() == ["A", "B", "C"]
        # v1 binds four times more weakly in C than in A; the three fits follow that order, which
        # a pooled fit could not report at all.
        kds = v1["kd"].to_list()
        assert kds == sorted(kds)
        assert kds[2] > kds[0] * 2

    def test_replicates_of_one_library_produce_one_value(self):
        outputs = run(_frame(IDENTICAL), params=LF, replicates=REPS)
        v1 = outputs["per_replicate"].filter(pl.col("clonotypeKey") == "v1")["kd"].to_list()
        assert v1[0] == pytest.approx(v1[1]) == pytest.approx(v1[2])

    def test_the_combined_value_is_the_mean_of_the_replicates_log10_kd(self):
        outputs = run(_frame(THREE), params=LF, replicates=REPS)
        per_replicate = outputs["per_replicate"].filter(pl.col("clonotypeKey") == "v1")
        expected = math.fsum(math.log10(v) for v in per_replicate["kd"].to_list()) / 3
        assert math.log10(_kd(outputs["per_clonotype"], "v1")) == pytest.approx(expected)

    def test_identical_replicates_report_a_zero_spread(self):
        outputs = run(_frame(IDENTICAL), params=LF, replicates=REPS)
        assert _kd(outputs["per_clonotype"], "v2", "kdLog10Sd") == pytest.approx(0.0, abs=1e-12)
        assert _kd(outputs["per_clonotype"], "v2", "kdLog10Se") == pytest.approx(0.0, abs=1e-12)

    def test_a_spread_across_replicates_becomes_the_reported_error(self):
        outputs = run(_frame(THREE), params=LF, replicates=REPS)
        row = outputs["per_clonotype"].filter(pl.col("clonotypeKey") == "v1").row(0, named=True)
        assert row["kdLog10Sd"] > 0.1
        assert row["kdLog10Se"] == pytest.approx(row["kdLog10Sd"] / math.sqrt(3))

    def test_the_counts_name_what_each_replicate_produced(self):
        design = {r: {"v1": 1e-8, "tight": 1e-14} for r in ("A", "B", "C")}
        outputs = run(_frame(design), params=LF, replicates=REPS)
        row = outputs["per_clonotype"].filter(pl.col("clonotypeKey") == "tight").row(0, named=True)
        assert row["replicateCount"] == 3
        assert row["replicateCensoredBelowCount"] == 3
        assert row["replicateMeasuredCount"] == 0
        assert row["kdOutOfRange"] is True

    def test_one_replicate_reproduces_a_single_titration(self):
        design = {"A": {"v1": 1e-8, "v2": 3e-9}}
        reads = _frame(design)
        combined = run(reads, params=LF, replicates=REPS)["per_clonotype"]
        plain = run(reads.drop("replicate"), params=LF)["per_clonotype"]
        for key in ("v1", "v2"):
            assert _kd(combined, key) == pytest.approx(_kd(plain, key), rel=1e-12)
            assert _kd(combined, key, "kdLog10Se") == pytest.approx(_kd(plain, key, "kdLog10Se"))

    def test_an_unusable_replicate_column_warns_and_runs_as_one_replicate(self, capsys):
        reads = _frame(THREE).with_columns(pl.lit(None, dtype=pl.Utf8).alias("replicate"))
        outputs = run(reads, params=LF, replicates=REPS)
        assert "WARN" in capsys.readouterr().out
        assert outputs["per_replicate"]["replicate"].unique().to_list() == ["1"]

    def test_the_signal_frames_keep_each_replicate_apart(self):
        outputs = run(_frame(THREE), params=LF, replicates=REPS)
        assert outputs["mean_bin"].columns == ["clonotypeKey", "replicate", "concentrationStr", "meanBin"]
        assert outputs["mean_bin"].filter(pl.col("clonotypeKey") == "v1")["replicate"].n_unique() == 3
        assert outputs["fitted_mean_bin"].columns == [
            "clonotypeKey",
            "replicate",
            "concentrationStr",
            "fittedMeanBin",
        ]

    def test_the_legacy_equation_also_fits_each_replicate(self):
        outputs = run(_frame(THREE), params=HILL, replicates=REPS)
        assert outputs["per_replicate"].height == 6
        assert _kd(outputs["per_clonotype"], "v1") is not None


class TestSortFractionAndDepth:
    def test_replicates_sharing_a_concentration_no_longer_break_the_sort_fraction_check(self):
        reads = _frame(THREE)
        with pytest.raises(InputValidationError, match="sum to 1.0"):
            run(reads.drop("replicate"), params=LF, sort_fraction_col="sort_fraction")
        outputs = run(reads, params=LF, sort_fraction_col="sort_fraction", replicates=REPS)
        assert outputs["per_clonotype"].height == 2

    def test_a_deeper_replicate_does_not_move_another_replicates_signal(self):
        design = {"A": {"v1": 1e-8}, "B": {"v1": 1e-8}}
        shallow = run(_frame({"A": design["A"]}), params=LF)["mean_bin"]
        together = run(_frame(design, B=8000), params=LF, replicates=REPS)["mean_bin"]
        a_only = together.filter(pl.col("replicate") == "A").drop("replicate")
        assert_frame_equal(a_only, shallow)


class TestFiltersEndToEnd:
    def _design(self):
        # "noisy" is a shallow, broad copy in replicate B: a wide fit error and a poor R².
        return {
            "A": {"v1": 1e-8},
            "B": {"v1": 1e-8},
        }

    def test_no_filter_keeps_every_replicate(self):
        outputs = run(_frame(self._design()), params=LF, replicates=REPS)
        assert outputs["per_clonotype"]["replicateFilteredCount"].to_list() == [0]

    def test_an_unreachable_fit_error_cutoff_filters_every_replicate(self):
        inputs = ReplicateInputs(column="replicate", filters=Filters(fit_error_cutoff=1e-9))
        outputs = run(_frame(self._design()), params=LF, replicates=inputs)
        row = outputs["per_clonotype"].row(0, named=True)
        assert row["replicateFilteredCount"] == 2
        assert row["affinityClass"] == "Failed"
        assert row["fitFailureReason"] == "all_replicates_filtered"
        assert row["kd"] is None

    def test_an_unreachable_r_squared_floor_filters_every_replicate(self):
        inputs = ReplicateInputs(column="replicate", filters=Filters(r2_floor=1.5))
        outputs = run(_frame(self._design()), params=LF, replicates=inputs)
        assert outputs["per_clonotype"]["fitFailureReason"].to_list() == ["all_replicates_filtered"]
        assert outputs["per_replicate"]["kdState"].unique().to_list() == ["Left out by a filter"]

    def test_a_reachable_cutoff_keeps_the_replicates(self):
        inputs = ReplicateInputs(column="replicate", filters=Filters(fit_error_cutoff=5.0, r2_floor=0.0))
        outputs = run(_frame(self._design()), params=LF, replicates=inputs)
        assert outputs["per_clonotype"]["replicateFilteredCount"].to_list() == [0]

    def test_the_mean_bin_fit_ignores_a_cutoff_it_cannot_judge(self):
        # The mean bin fit writes no per-replicate kdLog10Se. Reading that None as a rejection
        # rated every variant Failed on any run that merely had the cutoff set.
        inputs = ReplicateInputs(column="replicate", filters=Filters(fit_error_cutoff=1.0))
        outputs = run(_frame(self._design()), params=HILL, replicates=inputs)
        row = outputs["per_clonotype"].row(0, named=True)
        assert row["replicateFilteredCount"] == 0
        assert row["replicateCount"] == 2
        assert row["fitFailureReason"] != "all_replicates_filtered"


# Two clonotypes, each titrated in replicates A and B; "c2" only in A.
TWO_REPLICATES = pl.DataFrame(
    [
        row
        for nt, kd, rep in (
            ("c1", 1e-8, "A"),
            ("c2", 1e-7, "A"),
            ("c1", 1e-8, "B"),
        )
        for row in _reads(nt, kd, rep)
    ]
)


class TestReplicateCli:
    def test_the_cli_writes_the_per_replicate_frame(self, tmp_path):
        TWO_REPLICATES.write_parquet(tmp_path / "reads.parquet")
        (tmp_path / "params.json").write_text(json.dumps({"equation": "log_fluorescence"}))
        out = tmp_path / "out"
        out.mkdir()
        rc = main(
            [
                "--reads", str(tmp_path / "reads.parquet"),
                "--params", str(tmp_path / "params.json"),
                "--out-per-clonotype", str(out / "per_clonotype.tsv"),
                "--out-mean-bin", str(out / "mean_bin.tsv"),
                "--out-fitted-mean-bin", str(out / "fitted_mean_bin.tsv"),
                "--out-concentration-value", str(out / "concentration_value.tsv"),
                "--replicate-column", "replicate",
                "--replicate-error-model", "phillips",
                "--out-per-replicate", str(out / "per_replicate.tsv"),
            ]
        )  # fmt: skip
        assert rc == 0
        per_replicate = pl.read_csv(out / "per_replicate.tsv", separator="\t")
        assert per_replicate.height == 3
        per_clonotype = pl.read_csv(out / "per_clonotype.tsv", separator="\t")
        assert per_clonotype.filter(pl.col("clonotypeKey") == "c1")["replicateCount"].item() == 2

    def test_the_cli_refuses_a_censored_majority_below_a_half(self, tmp_path):
        with pytest.raises(SystemExit):
            main(
                [
                    "--reads", "r.parquet",
                    "--out-per-clonotype", "a", "--out-mean-bin", "b",
                    "--out-fitted-mean-bin", "c", "--out-concentration-value", "d",
                    "--replicate-column", "replicate", "--out-per-replicate", "e",
                    "--replicate-censored-majority", "0.2",
                ]
            )  # fmt: skip


class TestReviewClassesAcrossReplicates:
    DESIGN = {r: {"dim": 1e3, "v1": 1e-8} for r in ("A", "B")}

    def test_a_flat_dim_variant_in_every_replicate_is_above_range(self):
        per_clonotype = run(_frame(self.DESIGN), params=LF, replicates=REPS)["per_clonotype"]
        row = per_clonotype.filter(pl.col("clonotypeKey") == "dim").row(0, named=True)
        assert (row["affinityClass"], row["kdRangeSide"]) == ("Above range", "Above range")
        assert _kd(per_clonotype, "v1", "kdRangeSide") == "In range"

    def test_the_r2_gate_fails_each_replicate_inside_its_own_fit(self):
        gated = FitParams(equation="log_fluorescence", r2_gate_min=1.0)
        outputs = run(_frame(self.DESIGN), params=gated, replicates=REPS)
        per_replicate = outputs["per_replicate"].filter(pl.col("clonotypeKey") == "v1")
        assert per_replicate["kdState"].to_list() == ["Left out, no fit", "Left out, no fit"]
        assert per_replicate["fitFailureReason"].to_list() == ["low_r2", "low_r2"]
        assert _kd(outputs["per_clonotype"], "v1", "affinityClass") == "Failed"
        assert _kd(outputs["per_clonotype"], "v1", "fitFailureReason") == "all_replicates_failed"
