"""Sample replicates through the whole pipeline, on a synthetic library with known K_D.

Each replicate is a full titration of its own: its own samples, its own gate depths and its own
sort fractions. The expectations are the mean of the per-replicate fits, the counts, and — for a
run with no replicate column — that nothing moves at all.
"""

from __future__ import annotations

import json
import math
import statistics

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from constants import FitParams
from io_layer import InputValidationError
from main import main
from pipeline import run
from replicate_combine import COUNT_SCHEMA, ReplicateInputs
from synonymous import SynonymousInputs

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
    """With no replicate column asked for, every output is that of a run without replicates."""

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
                    "--out-per-replicate", "e",
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

    def test_identical_replicates_report_a_zero_spread_and_their_fit_error_as_the_se(self):
        outputs = run(_frame(IDENTICAL), params=LF, replicates=REPS)
        assert _kd(outputs["per_clonotype"], "v2", "kdLog10Sd") == pytest.approx(0.0, abs=1e-12)
        # The spread is zero, so the fit-error term sets the SE: RMS of the replicates' fit SEs / sqrt(3).
        fit_ses = outputs["per_replicate"].filter(pl.col("clonotypeKey") == "v2")["kdLog10Se"].to_list()
        assert len(fit_ses) == 3 and all(se > 0 for se in fit_ses)
        expected = math.sqrt(math.fsum(se * se for se in fit_ses) / 3) / math.sqrt(3)
        assert _kd(outputs["per_clonotype"], "v2", "kdLog10Se") == pytest.approx(expected, rel=1e-9)

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

    def test_each_replicate_reads_its_own_gate_range(self):
        # Replicate B was gated higher and narrower. Pooled with A, B's fit would be told its gates
        # span A's floor too, which moves its plateau threshold and its fit bounds.
        kds = {f"k{i}": kd for i, kd in enumerate([1e-7, 2e-7, 3e-7, 5e-7, 1e-6, 5e-6])}
        stretch = {"A": 1.0, "B": 0.75}
        offset = {"A": 0.0, "B": 0.5}
        reads = _frame({"A": kds, "B": kds}).with_columns(
            (
                (pl.col("bin") - 1) * GATE_WIDTH * pl.col("replicate").replace_strict(stretch, return_dtype=pl.Float64)
                + 1.0
                + pl.col("replicate").replace_strict(offset, return_dtype=pl.Float64)
            ).alias("bin_fluorescence")
        )
        combined = run(reads, params=LF, replicates=REPS)["per_replicate"]
        for replicate in ("A", "B"):
            alone = run(reads.filter(pl.col("replicate") == replicate).drop("replicate"), params=LF)["per_clonotype"]
            fitted = combined.filter(pl.col("replicate") == replicate)
            for key in kds:
                assert _kd(fitted, key) == pytest.approx(_kd(alone, key), rel=1e-9), (replicate, key)
                assert _kd(fitted, key, "affinityClass") == _kd(alone, key, "affinityClass"), (replicate, key)

    def test_the_legacy_equation_also_fits_each_replicate(self):
        outputs = run(_frame(THREE), params=HILL, replicates=REPS)
        assert outputs["per_replicate"].height == 6
        assert _kd(outputs["per_clonotype"], "v1") is not None


class TestDifferentConcentrationGrids:
    """Replicate A titrated without the top concentration, B without the bottom one."""

    @pytest.fixture(scope="class")
    def outputs(self):
        reads = _frame({"A": {"v1": 1e-8}, "B": {"v1": 2e-8}})
        reads = reads.filter(
            ~((pl.col("replicate") == "A") & (pl.col("concentration") == CONCS[-1]))
            & ~((pl.col("replicate") == "B") & (pl.col("concentration") == CONCS[0]))
        )
        return run(reads, params=LF, replicates=REPS)

    def test_every_replicate_curve_is_written_on_the_union_of_the_grids(self, outputs):
        fitted = outputs["fitted_mean_bin"]
        for replicate in ("A", "B"):
            concs = fitted.filter(pl.col("replicate") == replicate)["concentrationStr"].to_list()
            assert sorted(concs) == sorted(str(c) for c in CONCS)

    def test_each_replicate_fit_reads_only_its_own_points(self, outputs):
        alone = run(
            _frame({"A": {"v1": 1e-8}}).filter(pl.col("concentration") != CONCS[-1]).drop("replicate"), params=LF
        )
        a = outputs["per_replicate"].filter(pl.col("replicate") == "A")
        assert _kd(a, "v1") == pytest.approx(_kd(alone["per_clonotype"], "v1"), rel=1e-9)


class TestFitReadsWithReplicates:
    def test_a_replicate_gated_before_its_fit_adds_no_read(self):
        # v1 has every concentration in A; in B only two, under the 5-point gate.
        reads = _frame({"A": {"v1": 1e-8}, "B": {"v1": 1e-8}})
        reads = reads.filter((pl.col("replicate") == "A") | pl.col("concentration").is_in(CONCS[:2]))
        outputs = run(reads, params=LF, replicates=REPS, read_depth=True)
        b = outputs["per_replicate"].filter(pl.col("replicate") == "B").row(0, named=True)
        assert b["fitFailureReason"] == "insufficient_points"
        a_reads = reads.filter(pl.col("replicate") == "A")["reads"].sum()
        assert dict(outputs["fit_reads"].iter_rows())["v1"] == a_reads


class TestSortFractionAndDepth:
    def test_replicates_sharing_a_concentration_pass_the_sort_fraction_check(self):
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


class TestReplicateCounts:
    COUNTS = [
        "replicateCount",
        "replicateMeasuredCount",
        "replicateCensoredBelowCount",
        "replicateCensoredAboveCount",
        "replicateFailedCount",
    ]

    def test_every_replicate_enters_the_mean(self):
        outputs = run(_frame({"A": {"v1": 1e-8}, "B": {"v1": 1e-8}}), params=LF, replicates=REPS)
        row = outputs["per_clonotype"].row(0, named=True)
        assert [row[c] for c in self.COUNTS] == [2, 2, 0, 0, 0]

    def test_both_combined_frames_carry_exactly_these_counts(self):
        assert list(COUNT_SCHEMA) == self.COUNTS
        outputs = run(UNBALANCED, params=LF, synonymous=SynonymousInputs(nt_to_aa=UNBALANCED_LINKS), replicates=REPS)
        for name in ("per_clonotype", "per_variant"):
            counts = {c for c in outputs[name].columns if c.startswith("replicate") and c.endswith("Count")}
            assert counts == set(self.COUNTS), name


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
                "--out-per-replicate", str(out / "per_replicate.tsv"),
            ]
        )  # fmt: skip
        assert rc == 0
        per_replicate = pl.read_csv(out / "per_replicate.tsv", separator="\t")
        assert per_replicate.height == 3
        per_clonotype = pl.read_csv(out / "per_clonotype.tsv", separator="\t")
        assert per_clonotype.filter(pl.col("clonotypeKey") == "c1")["replicateCount"].item() == 2

    def test_the_cli_refuses_the_nucleotide_per_replicate_frame_without_the_linker(self):
        with pytest.raises(SystemExit):
            main(
                [
                    "--reads", "r.parquet",
                    "--out-per-clonotype", "a", "--out-mean-bin", "b",
                    "--out-fitted-mean-bin", "c", "--out-concentration-value", "d",
                    "--replicate-column", "replicate", "--out-per-replicate", "e",
                    "--out-per-replicate-nt", "f",
                ]
            )  # fmt: skip

    @pytest.mark.parametrize(
        "option, value",
        [
            ("--replicate-error-model", "sd_over_sqrt_n"),
            ("--replicate-fit-error-cutoff", "1.0"),
            ("--replicate-r2-floor", "0.8"),
            ("--replicate-censored-majority", "0.5"),
        ],
    )
    def test_the_cli_has_no_replicate_option_beyond_the_column(self, option, value):
        with pytest.raises(SystemExit):
            main(
                [
                    "--reads", "r.parquet",
                    "--out-per-clonotype", "a", "--out-mean-bin", "b",
                    "--out-fitted-mean-bin", "c", "--out-concentration-value", "d",
                    "--replicate-column", "replicate", "--out-per-replicate", "e",
                    option, value,
                ]
            )  # fmt: skip


class TestReviewClassesAcrossReplicates:
    DESIGN = {r: {"dim": 1e3, "v1": 1e-8} for r in ("A", "B")}

    def test_a_flat_dim_variant_in_every_replicate_is_above_range(self):
        per_clonotype = run(_frame(self.DESIGN), params=LF, replicates=REPS)["per_clonotype"]
        row = per_clonotype.filter(pl.col("clonotypeKey") == "dim").row(0, named=True)
        assert (row["affinityClass"], row["kdRangeSide"]) == ("Above range", "Above range")
        assert _kd(per_clonotype, "v1", "kdRangeSide") == "In range"

    def test_the_combined_r2_reads_only_in_range_replicates(self):
        # One replicate in range, one pinned above: the combined R² is the in-range replicate's own,
        # and a variant pinned in every replicate has none.
        design = {"A": {"mix": 1e-8, "dim": 1e3}, "B": {"mix": 1e3, "dim": 1e3}}
        outputs = run(_frame(design), params=LF, replicates=REPS)
        per_replicate = outputs["per_replicate"].filter(pl.col("clonotypeKey") == "mix").sort("replicate")
        in_range_r2 = per_replicate["r2"].to_list()[0]
        assert per_replicate["r2"].to_list()[1] is None
        assert _kd(outputs["per_clonotype"], "mix", "r2") == pytest.approx(in_range_r2)
        assert _kd(outputs["per_clonotype"], "dim", "r2") is None

    def test_the_r2_gate_fails_each_replicate_inside_its_own_fit(self):
        gated = FitParams(equation="log_fluorescence", r2_gate_min=1.0)
        outputs = run(_frame(self.DESIGN), params=gated, replicates=REPS)
        per_replicate = outputs["per_replicate"].filter(pl.col("clonotypeKey") == "v1")
        assert per_replicate["affinityClass"].to_list() == ["Failed", "Failed"]
        assert per_replicate["fitFailureReason"].to_list() == ["low_r2", "low_r2"]
        assert _kd(outputs["per_clonotype"], "v1", "affinityClass") == "Failed"
        assert _kd(outputs["per_clonotype"], "v1", "fitFailureReason") == "all_replicates_failed"


# One amino-acid variant with two synonymous variants in replicate A and one in replicate B. The
# synonymous variants' K_D differ, so Starr's order (combine them, then the replicates) and a flat mean
# over all three fits do not agree.
UNBALANCED = pl.DataFrame(
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
UNBALANCED_LINKS = pl.DataFrame({"ntKey": ["c1", "c2"], "aaKey": ["X", "X"]})


class TestStarrOrdering:
    @pytest.fixture(scope="class")
    def outputs(self):
        return run(UNBALANCED, params=LF, synonymous=SynonymousInputs(nt_to_aa=UNBALANCED_LINKS), replicates=REPS)

    def test_the_synonymous_variants_combine_inside_each_replicate_before_the_replicates_combine(self, outputs):
        # The nt fits do not depend on the linker, so a run without it gives each nt variant's own value.
        nt = run(UNBALANCED, params=LF, replicates=REPS)["per_replicate"]
        log10 = {(r["clonotypeKey"], r["replicate"]): math.log10(r["kd"]) for r in nt.iter_rows(named=True)}
        starr = ((log10[("c1", "A")] + log10[("c2", "A")]) / 2 + log10[("c1", "B")]) / 2
        flat = (log10[("c1", "A")] + log10[("c2", "A")] + log10[("c1", "B")]) / 3
        got = math.log10(outputs["per_variant"].filter(pl.col("aaKey") == "X")["kd"].item())
        assert got == pytest.approx(starr)
        assert got != pytest.approx(flat)

    def test_a_stop_variant_fails_for_its_stop_codon_not_for_its_replicates(self):
        productive = pl.DataFrame({"aaKey": ["X"], "productive": ["false"]})
        outputs = run(
            UNBALANCED,
            params=LF,
            synonymous=SynonymousInputs(nt_to_aa=UNBALANCED_LINKS, aa_productive=productive),
            replicates=REPS,
        )
        x = outputs["per_variant"].filter(pl.col("aaKey") == "X").row(0, named=True)
        assert (x["kd"], x["affinityClass"], x["fitFailureReason"]) == (None, "Failed", "not_productive")
        per_replicate = outputs["per_replicate"].filter(pl.col("aaKey") == "X")
        assert set(per_replicate["fitFailureReason"].to_list()) == {"not_productive"}

    def test_the_per_replicate_frame_is_on_the_amino_acid_grain(self, outputs):
        frame = outputs["per_replicate"]
        assert frame.columns[:2] == ["aaKey", "replicate"]
        assert sorted(frame["replicate"].to_list()) == ["A", "B"]
        # Replicate A's reads are both synonymous variants'.
        reads = dict(zip(frame["replicate"].to_list(), frame["replicateReads"].to_list()))
        assert reads["A"] > reads["B"]

    def test_the_per_replicate_frame_carries_the_synonymous_variants_inside_each_replicate(self, outputs):
        nt = run(UNBALANCED, params=LF, replicates=REPS)["per_replicate"]
        log10 = {(r["clonotypeKey"], r["replicate"]): math.log10(r["kd"]) for r in nt.iter_rows(named=True)}
        rows = {r["replicate"]: r for r in outputs["per_replicate"].iter_rows(named=True)}
        assert rows["A"]["kdLog10SynonymSd"] == pytest.approx(
            statistics.stdev([log10[("c1", "A")], log10[("c2", "A")]])
        )
        assert (rows["A"]["ntVariants"], rows["A"]["ntMeasured"]) == (2, 2)
        # One synonymous variant has no spread.
        assert rows["B"]["kdLog10SynonymSd"] is None
        assert (rows["B"]["ntVariants"], rows["B"]["ntMeasured"]) == (1, 1)

    def test_the_nucleotide_per_replicate_frame_holds_every_nt_variant_in_every_replicate(self, outputs):
        nt = run(UNBALANCED, params=LF, replicates=REPS)["per_replicate"]
        frame = outputs["per_replicate_nt"]
        assert frame.columns[:2] == ["clonotypeKey", "replicate"]
        assert (
            frame.sort("clonotypeKey", "replicate")
            .select("clonotypeKey", "replicate", "kd")
            .equals(nt.sort("clonotypeKey", "replicate").select("clonotypeKey", "replicate", "kd"))
        )

    def test_each_replicate_s_synonymous_se_is_the_fit_se_the_replicates_combine(self, outputs):
        per_replicate = {r["replicate"]: r for r in outputs["per_replicate"].iter_rows(named=True)}
        values = [math.log10(per_replicate[r]["kd"]) for r in ("A", "B")]
        ses = [per_replicate[r]["kdLog10Se"] for r in ("A", "B")]
        spread = statistics.stdev(values) / math.sqrt(2)
        fit_term = math.sqrt((ses[0] ** 2 + ses[1] ** 2) / 2) / math.sqrt(2)
        combined = outputs["per_variant"].filter(pl.col("aaKey") == "X")["kdLog10Se"].item()
        assert combined == pytest.approx(max(spread, fit_term), rel=1e-9)

    def test_the_nucleotide_counts_add_up_over_the_replicates(self, outputs):
        row = outputs["per_variant"].row(0, named=True)
        assert row["ntVariants"] == 3
        assert row["ntMeasured"] == 3
        assert row["replicateCount"] == 2

    def test_the_summary_names_the_replicates_and_the_error_model(self, outputs):
        summary = outputs["synonymous_summary"]
        assert summary["replicates"] == 2
        assert summary["errorModel"] == "sd_over_sqrt_n"
        assert summary["aaVariants"] == 1
        assert summary["ntVariants"] == 2


# Two replicates, each with two synonymous variants of close K_D.
AGREEING = pl.DataFrame(
    [row for rep in ("A", "B") for nt, kd in (("c1", 1e-8), ("c2", 1.1e-8)) for row in _reads(nt, kd, rep)]
)
AGREEING_LINKS = pl.DataFrame({"ntKey": ["c1", "c2"], "aaKey": ["X", "X"]})


class TestGoodMajorityUnderTheMeanBinFit:
    """Under "hill_gate_index" each replicate's amino-acid value carries the spread of its synonymous
    variants as its SE. That SE is no fit SE, so the replicates' mean still needs a Good majority."""

    # No R² reaches 1, so every in-range fit is Partial.
    PARTIAL_FITS = FitParams(equation="hill_gate_index", r2_threshold_good=1.0)

    def _outputs(self):
        return run(
            AGREEING, params=self.PARTIAL_FITS, synonymous=SynonymousInputs(nt_to_aa=AGREEING_LINKS), replicates=REPS
        )

    def test_agreeing_partial_fits_combine_into_a_partial_amino_acid_value(self):
        outputs = self._outputs()
        assert set(outputs["per_replicate_nt"]["affinityClass"].to_list()) == {"Partial"}
        per_replicate = outputs["per_replicate"]
        assert per_replicate["affinityClass"].to_list() == ["Partial", "Partial"]
        assert all(se is not None for se in per_replicate["kdLog10Se"].to_list())
        row = outputs["per_variant"].row(0, named=True)
        assert row["kdLog10Se"] <= 0.5
        assert row["affinityClass"] == "Partial"
