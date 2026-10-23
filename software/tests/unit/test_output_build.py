"""Behavioral tests for output_build.py (R14/R14b)."""

from __future__ import annotations

import polars as pl
import pytest
from hypothesis import given
from hypothesis import strategies as st

from output_build import (
    PER_CLONOTYPE_SCHEMA,
    add_diagnostic_plot_columns,
    build_concentration_value_frame,
    build_mean_bin_frame,
    flag_kd_out_of_range,
)


def _per_clonotype(*kds: float | None) -> pl.DataFrame:
    """Minimal per-clonotype frame, one row per Kd; a null Kd is a Failed fit with Hill coefficient 1.0."""
    return pl.DataFrame(
        [
            {
                "clonotypeKey": f"C{i}",
                "kd": kd,
                "hillCoefficient": 1.0,
                "r2": 0.9,
                "affinityClass": "Good" if kd is not None else "Failed",
                "fitFailureReason": None if kd is not None else "low_r2",
                "kdOutOfRange": None,
            }
            for i, kd in enumerate(kds)
        ],
        schema=PER_CLONOTYPE_SCHEMA,
    )


positive_kd = st.floats(1e-15, 1e15)


class TestFlagKdOutOfRange:
    # R14b: outside [min, max] → kdOutOfRange = True. Closed interval (boundary = in-range).
    @given(kds=st.lists(positive_kd, min_size=1, max_size=20), lo=positive_kd, span=st.floats(1.0, 1e6))
    def test_flag_is_true_exactly_outside_the_range(self, kds, lo, span):
        hi = lo * span
        out = flag_kd_out_of_range(_per_clonotype(*kds), min_concentration=lo, max_concentration=hi)
        assert out["kdOutOfRange"].to_list() == [not lo <= kd <= hi for kd in kds]

    # The closed interval is easy to break with < vs <=, and generated values almost never land on it.
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

    @given(kds=st.lists(st.none() | positive_kd, min_size=1, max_size=20), max_c=positive_kd)
    def test_plot_positions_keep_fitted_values_and_park_null_kd_one_decade_right(self, kds, max_c):
        out = add_diagnostic_plot_columns(_per_clonotype(*kds), max_concentration=max_c)
        assert out["kdPlotPosition"].to_list() == [max_c * 10.0 if kd is None else kd for kd in kds]
        assert out["hillPlotPosition"].to_list() == [1.0] * len(kds)
        assert out.schema["kdPlotPosition"] == pl.Float64
        assert out.schema["hillPlotPosition"] == pl.Float64

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
