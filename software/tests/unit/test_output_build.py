"""Behavioral tests for output_build.py (R14/R14b)."""

from __future__ import annotations

import polars as pl
import pytest

from output_build import (
    PER_CLONOTYPE_SCHEMA,
    add_diagnostic_plot_columns,
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
    # R17: Failed clonotypes with null K_D must appear at a sentinel x right of the fitted range
    # so the Affinity-vs-Fit scatter doesn't silently drop them on a log axis.

    def test_null_kd_maps_to_decade_right_of_max(self):
        frame = add_diagnostic_plot_columns(_per_clonotype(None), max_concentration=100.0)
        assert frame["kdPlotPosition"][0] == 1000.0

    def test_finite_kd_passes_through_unchanged(self):
        frame = add_diagnostic_plot_columns(_per_clonotype(5.0), max_concentration=100.0)
        assert frame["kdPlotPosition"][0] == 5.0

    def test_null_hill_maps_to_one(self):
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
        assert out["hillPlotPosition"][0] == 1.0

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
    def test_c0_excluded(self):
        # TiteSeq assay range is sub-µM; use 100 nM so the concentration stays
        # well below the attomolar-encoding ceiling (~9.22 M) enforced by the
        # output-build Int64 cast and the R2 validation guard.
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

    def test_preserves_canonical_string(self):
        signal = pl.DataFrame(
            [
                {"clonotypeKey": "A", "concentrationStr": "1e-6", "concentration": 1e-6, "signal": 2.0},
            ]
        )
        out = build_mean_bin_frame(signal)
        assert out["concentrationStr"][0] == "1e-6"

    def test_equivalent_strings_share_concentrationAM(self):
        # R14: two different canonical strings that parse to the same float MUST yield
        # the same attomolar integer key, so the numeric sibling axis cannot drift away
        # from the canonical string axis for equal concentrations.
        signal = pl.DataFrame(
            [
                {"clonotypeKey": "A", "concentrationStr": "1e-7", "concentration": 1e-7, "signal": 2.0},
                {"clonotypeKey": "B", "concentrationStr": "0.0000001", "concentration": 1e-7, "signal": 3.0},
            ]
        )
        out = build_mean_bin_frame(signal)
        assert out["concentrationAM"].to_list() == [100_000_000_000, 100_000_000_000]
        # Canonical strings are preserved distinct even though the numeric axis matches.
        assert out["concentrationStr"].to_list() == ["1e-7", "0.0000001"]
