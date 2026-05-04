"""Behavioral tests for output_build.py (R14/R14b)."""

from __future__ import annotations

import polars as pl
import pytest

from output_build import (
    PER_CLONOTYPE_SCHEMA,
    add_diagnostic_plot_columns,
    build_concentration_value_frame,
    build_mean_bin_frame,
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


class TestAddDiagnosticPlotColumns:
    # R17: Failed clonotypes with null Kd must appear at a sentinel x right of the fitted range
    # so the Affinity-vs-Fit scatter doesn't silently drop them on a log axis.

    def test_null_kd_maps_to_decade_right_of_max(self):
        frame = add_diagnostic_plot_columns(_per_clonotype(None), max_concentration=100.0)
        assert frame["kdPlotPosition"][0] == 1000.0

    def test_finite_kd_passes_through_unchanged(self):
        frame = add_diagnostic_plot_columns(_per_clonotype(5.0), max_concentration=100.0)
        assert frame["kdPlotPosition"][0] == 5.0

    def test_null_hill_maps_to_minus_one(self):
        # Sentinel is -1.0: non-physical (Hill is strictly positive), so Failed rows
        # pool visibly below the n>0 cluster instead of mixing with well-fitted n≈1.
        frame = pl.DataFrame(
            [
                {
                    "clonotypeKey": "A",
                    "kd": None,
                    "hillCoefficient": None,
                    "r2": None,
                    "affinityClass": "Failed",
                    "fitFailureReason": "convergence_failure",
                    "kdOutOfRange": None,
                }
            ],
            schema=PER_CLONOTYPE_SCHEMA,
        )
        out = add_diagnostic_plot_columns(frame, max_concentration=100.0)
        assert out["hillPlotPosition"][0] == -1.0

    def test_finite_hill_passes_through_unchanged(self):
        out = add_diagnostic_plot_columns(_per_clonotype(5.0), max_concentration=100.0)
        # _per_clonotype sets hillCoefficient = 1.0 for non-null kd rows
        assert out["hillPlotPosition"][0] == 1.0

    def test_output_columns_are_float64_and_never_null(self):
        frame = pl.concat([_per_clonotype(None), _per_clonotype(5.0)])
        out = add_diagnostic_plot_columns(frame, max_concentration=100.0)
        assert out.schema["kdPlotPosition"] == pl.Float64
        assert out.schema["hillPlotPosition"] == pl.Float64
        assert out["kdPlotPosition"].null_count() == 0
        assert out["hillPlotPosition"].null_count() == 0


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
        signal = pl.DataFrame(
            [{"clonotypeKey": "A", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 2.0}]
        )
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
