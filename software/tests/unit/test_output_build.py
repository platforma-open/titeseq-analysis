"""Behavioral tests for output_build.py (R14/R14b)."""

from __future__ import annotations

import math

import polars as pl
import pytest

from constants import COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL
from output_build import (
    COL_AA_KEY,
    FITTED_MEAN_BIN_SCHEMA,
    PER_CLONOTYPE_SCHEMA,
    add_log10_kd,
    build_aa_curve_frame,
    build_concentration_value_frame,
    build_mean_bin_frame,
    build_variant_curve_frame,
    curve_variants,
    flag_kd_out_of_range,
)


def _per_clonotype(kd):
    """Minimal per-clonotype frame with one row for flag_kd_out_of_range tests."""
    return pl.DataFrame(
        [
            {
                "clonotypeKey": "A",
                "kd": kd,
                "hillCoefficient": 1.0,
                "r2": 0.9,
                "affinityClass": "Good" if kd is not None else "Failed",
                "fitFailureReason": None if kd is not None else "low_r2",
                "kdOutOfRange": None,
            }
        ],
        schema=PER_CLONOTYPE_SCHEMA,
    )


class TestFlagKdOutOfRange:
    # R14b: outside [min, max] → kdOutOfRange = True. Closed interval (boundary = in-range).
    def test_in_range_false(self):
        out = flag_kd_out_of_range(_per_clonotype(1.0), min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is False

    def test_above_max_flag_true(self):
        out = flag_kd_out_of_range(_per_clonotype(1000.0), min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is True

    def test_below_min_flag_true(self):
        out = flag_kd_out_of_range(_per_clonotype(0.01), min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is True

    @pytest.mark.parametrize("kd_value", [0.1, 100.0])
    def test_exactly_at_boundary_in_range(self, kd_value):
        out = flag_kd_out_of_range(_per_clonotype(kd_value), min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is False

    def test_failed_fit_kd_null_stays_null(self):
        out = flag_kd_out_of_range(_per_clonotype(None), min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is None


class TestMeanBinFrame:
    # R14 mandates the c=0 sample is dropped from output (log(0) = -inf breaks the graph).
    def test_c0_excluded(self):
        signal = pl.DataFrame(
            [
                {"clonotypeKey": "A", "concentrationStr": "0", "concentration": 0.0, "signal": 1.2},
                {"clonotypeKey": "A", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 2.5},
            ]
        )
        out = build_mean_bin_frame(signal)
        assert out.height == 1
        assert out["concentrationStr"][0] == "1e-7"
        assert out["meanBin"][0] == 2.5

    # Spec-realigned schema — clonotypeKey + concentrationStr (canonical String axis
    # key) + meanBin. No more concentrationAM, no Float concentration column.
    def test_output_columns_match_spec(self):
        signal = pl.DataFrame([{"clonotypeKey": "A", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 2.0}])
        out = build_mean_bin_frame(signal)
        assert out.columns == ["clonotypeKey", "concentrationStr", "meanBin"]

    # R14: the canonical string flows through to the output unchanged. The Tengo
    # workflow wraps it directly as a String axis on the meanBin PColumn; if it
    # drifts here, the axis keys mismatch downstream.
    @pytest.mark.parametrize(
        "conc_str",
        ["1e-12", "1e-10", "1e-9", "2.5e-9", "5e-8", "1e-7", "0.0000001", "1e-6"],
    )
    def test_canonical_string_preserved_byte_for_byte(self, conc_str):
        signal = pl.DataFrame(
            [{"clonotypeKey": "A", "concentrationStr": conc_str, "concentration": float(conc_str), "signal": 1.0}]
        )
        out = build_mean_bin_frame(signal)
        assert out["concentrationStr"][0] == conc_str

    # Regression: upstream `polars.group_by` in normalize() does not guarantee
    # row order, so the output must sort itself to keep the TSV byte-stable
    # across runs. Without the sort this drove CIDConflictError downstream.
    def test_output_sorted_by_clonotype_then_concentration(self):
        # Deliberately shuffled input — clonotypes interleaved, concentrations
        # not monotonic per clonotype.
        signal = pl.DataFrame(
            [
                {"clonotypeKey": "B", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 2.0},
                {"clonotypeKey": "A", "concentrationStr": "1e-9", "concentration": 1e-9, "signal": 1.0},
                {"clonotypeKey": "B", "concentrationStr": "1e-9", "concentration": 1e-9, "signal": 3.0},
                {"clonotypeKey": "A", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 4.0},
            ]
        )
        out = build_mean_bin_frame(signal)
        keys = list(zip(out["clonotypeKey"].to_list(), out["concentrationStr"].to_list()))
        assert keys == sorted(keys)


class TestConcentrationValueFrame:
    # The concentrationValue sidecar PColumn provides the numeric source for the
    # Titration Curves X-axis. Axes [concentration:String]; valueType Double; one
    # row per unique non-zero concentration.
    def test_dedup_per_concentration(self):
        signal = pl.DataFrame(
            [
                {"clonotypeKey": "A", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 2.0},
                {"clonotypeKey": "B", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 3.0},
                {"clonotypeKey": "A", "concentrationStr": "1e-6", "concentration": 1e-6, "signal": 4.0},
            ]
        )
        out = build_concentration_value_frame(signal)
        assert out.columns == ["concentrationStr", "concentration"]
        assert sorted(out["concentrationStr"].to_list()) == ["1e-6", "1e-7"]
        for row in out.iter_rows(named=True):
            assert row["concentration"] == pytest.approx(float(row["concentrationStr"]))

    def test_c0_excluded(self):
        # Same rationale as TestMeanBinFrame.test_c0_excluded — log(0) on the
        # log-scale X-axis is undefined.
        signal = pl.DataFrame(
            [
                {"clonotypeKey": "A", "concentrationStr": "0", "concentration": 0.0, "signal": 1.2},
                {"clonotypeKey": "A", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 2.5},
            ]
        )
        out = build_concentration_value_frame(signal)
        assert "0" not in out["concentrationStr"].to_list()


class TestNegLog10Kd:
    def test_it_is_minus_log10_of_the_reported_kd(self):
        frame = pl.DataFrame({"kd": [1e-9, 2.5e-8]})
        got = add_log10_kd(frame)["kdNegLog10"].to_list()
        assert got == pytest.approx([9.0, -math.log10(2.5e-8)])

    def test_a_tighter_binder_reads_larger(self):
        got = add_log10_kd(pl.DataFrame({"kd": [1e-6, 1e-10]}))["kdNegLog10"].to_list()
        assert got[1] > got[0]

    def test_a_clonotype_with_no_kd_has_none(self):
        frame = pl.DataFrame({"kd": [None]}, schema={"kd": pl.Float64})
        assert add_log10_kd(frame)["kdNegLog10"].to_list() == [None]

    def test_a_non_positive_kd_has_none_rather_than_minus_infinity(self):
        frame = pl.DataFrame({"kd": [0.0]})
        assert add_log10_kd(frame)["kdNegLog10"].to_list() == [None]


class TestLog10Kd:
    def test_it_is_log10_of_the_reported_kd(self):
        frame = pl.DataFrame({"kd": [1e-9, 2.5e-8]})
        got = add_log10_kd(frame)["kdLog10"].to_list()
        assert got == pytest.approx([-9.0, math.log10(2.5e-8)])

    def test_it_is_the_negation_of_minus_log10(self):
        out = add_log10_kd(pl.DataFrame({"kd": [3e-7, 1e-11]}))
        assert out["kdLog10"].to_list() == pytest.approx([-v for v in out["kdNegLog10"].to_list()])

    def test_no_kd_and_non_positive_kd_have_none(self):
        frame = pl.DataFrame({"kd": [None, 0.0]}, schema={"kd": pl.Float64})
        assert add_log10_kd(frame)["kdLog10"].to_list() == [None, None]


def _aa_curve_inputs():
    rows = [(nt, rep, c) for nt in ("n1", "n2", "n3") for rep in ("A", "B") for c in ("1e-9", "1e-8")]
    fitted_value = {"n1": 1.0, "n2": 3.0, "n3": 10.0}
    fitted = pl.DataFrame(
        {
            COL_CLONOTYPE: [r[0] for r in rows],
            "replicate": [r[1] for r in rows],
            COL_CONC_STR: [r[2] for r in rows],
            "fittedMeanBin": [fitted_value[r[0]] + (1.0 if r[1] == "B" else 0.0) for r in rows],
        }
    ).filter(~((pl.col(COL_CLONOTYPE) == "n2") & (pl.col("replicate") == "B")))
    link = pl.DataFrame({COL_CLONOTYPE: ["n1", "n2", "n3", "n9"], COL_AA_KEY: ["a1", "a1", "a2", "a2"]})
    conc = pl.DataFrame({COL_CONC_STR: ["1e-9", "1e-8"], COL_CONC_VAL: [1e-9, 1e-8]})
    return fitted, link, conc


class TestAaCurveFrame:
    def test_replicate_curve_is_the_mean_of_the_nt_curves(self):
        out = build_aa_curve_frame(*_aa_curve_inputs())
        a1 = out.filter((pl.col("aaKey") == "a1") & (pl.col(COL_CONC_STR) == "1e-9")).sort("replicate")
        assert a1["aaFittedMeanBin"].to_list() == [2.0, 2.0]  # A: mean(1, 3); B: n1 alone, 1 + 1
        # The mean over replicates repeats on each replicate's row.
        assert a1["aaCombinedFittedMeanBin"].to_list() == [2.0, 2.0]
        assert a1["aaSignalConcentration"].to_list() == [1e-9, 1e-9]

    def test_one_row_per_variant_replicate_and_concentration(self):
        out = build_aa_curve_frame(*_aa_curve_inputs())
        assert out.height == 8
        assert out.select("aaKey", "replicate", COL_CONC_STR).is_unique().all()

    def test_nothing_fitted_gives_no_rows(self):
        _, link, conc = _aa_curve_inputs()
        empty = pl.DataFrame({k: [] for k in FITTED_MEAN_BIN_SCHEMA}, schema=FITTED_MEAN_BIN_SCHEMA)
        assert build_aa_curve_frame(empty, link, conc).height == 0


class TestCurveVariants:
    def test_keeps_the_nt_variants_with_a_curve(self):
        fitted, link, _ = _aa_curve_inputs()
        out = curve_variants(link, fitted)
        assert out.rows() == [("n1", "a1"), ("n2", "a1"), ("n3", "a2")]


def _aa_per_variant(se: float | None = 0.2):
    return pl.DataFrame(
        {COL_AA_KEY: ["a1", "a2"], "kd": [3e-9, None], "kdLog10Se": [se, None]},
        schema={COL_AA_KEY: pl.Utf8, "kd": pl.Float64, "kdLog10Se": pl.Float64},
    )


class TestAaKdMarks:
    def _marks(self, per_variant=None):
        out = build_aa_curve_frame(*_aa_curve_inputs(), _aa_per_variant() if per_variant is None else per_variant)
        return out, out.filter(pl.col(COL_CONC_STR).str.starts_with("kd-"))

    def _bar_xs(self, per_variant=None):
        _, marks = self._marks(per_variant)
        bar = marks.filter(pl.col("aaKdSeRange").is_not_null() & (pl.col("replicate") == "A")).sort(COL_CONC_STR)
        return bar["aaSignalConcentration"].to_list()

    def test_the_line_stands_at_kd_over_the_height_of_the_curves(self):
        _, marks = self._marks()
        line = marks.filter(pl.col("aaKdLine").is_not_null() & (pl.col("replicate") == "A")).sort(COL_CONC_STR)
        assert line["aaSignalConcentration"].to_list() == pytest.approx([3e-9, 3e-9])
        assert line["aaKdLine"].to_list() == [2.0, 2.0]  # a1's replicate curves are both flat at 2

    def test_the_bar_runs_one_se_either_side_in_rising_x(self):
        xs = self._bar_xs()
        assert xs == sorted(xs) and len(set(xs)) == 6
        assert xs[0] == pytest.approx(3e-9 / 10**0.2) and xs[-1] == pytest.approx(3e-9 * 10**0.2)

    def test_the_bar_stops_at_the_concentrations_and_loses_the_cap_there(self):
        # One SE reaches 3e-10 and 3e-8; the concentrations run from 1e-9 to 1e-8.
        xs = self._bar_xs(_aa_per_variant(se=1.0))
        assert xs == sorted(xs) and len(xs) == 2
        assert xs[0] == pytest.approx(1e-9) and xs[-1] == pytest.approx(1e-8)

    def test_a_bar_cut_on_one_side_keeps_the_other_cap(self):
        # Half a decade reaches 9.5e-10, below the lowest concentration, and 9.5e-9, inside them.
        xs = self._bar_xs(_aa_per_variant(se=0.5))
        assert xs == sorted(xs) and len(xs) == 4
        assert xs[0] == pytest.approx(1e-9) and xs[-1] == pytest.approx(3e-9 * 10**0.5)

    def test_a_kd_outside_the_concentrations_keeps_its_side_of_the_bar(self):
        per_variant = _aa_per_variant(se=0.5).with_columns(
            pl.when(pl.col(COL_AA_KEY) == "a1").then(1e-7).otherwise(pl.col("kd")).alias("kd")
        )
        xs = self._bar_xs(per_variant)
        # Kd,app is above every concentration: the bar runs from one SE below it, capped, to Kd,app.
        assert xs == sorted(xs) and len(xs) == 4
        assert xs[0] == pytest.approx(1e-7 / 10**0.5) and xs[-1] == pytest.approx(1e-7)

    def test_marks_sit_on_every_replicate_and_only_on_variants_with_a_kd(self):
        _, marks = self._marks()
        assert marks["aaKey"].unique().to_list() == ["a1"]
        assert marks.group_by("replicate").len().sort("replicate")["len"].to_list() == [8, 8]

    def test_every_row_has_a_height_and_marks_carry_no_curve(self):
        out, marks = self._marks()
        assert out["aaCurveExtent"].null_count() == 0
        assert marks["aaFittedMeanBin"].null_count() == marks.height
        curves = out.filter(~pl.col(COL_CONC_STR).str.starts_with("kd-"))
        assert curves["aaCurveExtent"].to_list() == curves["aaCombinedFittedMeanBin"].to_list()

    def test_an_se_too_large_to_place_gives_the_line_alone(self):
        _, marks = self._marks(_aa_per_variant(se=1e6))
        assert marks["aaKdSeRange"].null_count() == marks.height

    def test_no_se_gives_the_line_alone(self):
        _, marks = self._marks(_aa_per_variant(se=None))
        assert marks["aaKdSeRange"].null_count() == marks.height
        assert marks.height == 4


class TestVariantCurveFrame:
    def _kd(self):
        return pl.DataFrame({COL_CLONOTYPE: ["n1"], "kd": [3e-9], "kdLog10Se": [0.2]})

    def test_with_replicates_the_variant_has_its_mean_curve_and_marks(self):
        fitted, _, conc = _aa_curve_inputs()
        out = build_variant_curve_frame(fitted, conc, self._kd())
        n1 = out.filter((pl.col(COL_CLONOTYPE) == "n1") & (pl.col(COL_CONC_STR) == "1e-9")).sort("replicate")
        assert n1["curveCombinedFittedMeanBin"].to_list() == [1.5, 1.5]
        assert n1["signalConcentration"].to_list() == [1e-9, 1e-9]
        marks = out.filter(pl.col(COL_CONC_STR).str.starts_with("kd-"))
        assert marks[COL_CLONOTYPE].unique().to_list() == ["n1"]
        assert marks.group_by("replicate").len()["len"].to_list() == [8, 8]

    def test_the_fitted_curve_is_not_repeated(self):
        fitted, _, conc = _aa_curve_inputs()
        out = build_variant_curve_frame(fitted, conc, self._kd())
        assert "curveFittedMeanBin" not in out.columns

    def test_without_replicates_there_is_no_replicate_column(self):
        fitted, _, conc = _aa_curve_inputs()
        one = fitted.filter(pl.col("replicate") == "A").drop("replicate")
        out = build_variant_curve_frame(one, conc, self._kd())
        assert "replicate" not in out.columns
        assert out.select(COL_CLONOTYPE, COL_CONC_STR).is_unique().all()
        line = out.filter(pl.col("kdLine").is_not_null()).sort(COL_CONC_STR)
        assert line["kdLine"].to_list() == [1.0, 1.0]
        assert out["curveExtent"].null_count() == 0
