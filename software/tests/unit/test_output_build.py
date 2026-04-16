"""Behavioral tests for output_build.py (R13/R14/R14b)."""

from __future__ import annotations

import polars as pl
import pytest

from output_build import (
    build_fitted_mean_bin_frame,
    build_mean_bin_frame,
    build_per_clonotype_frame,
    flag_kd_out_of_range,
)


class TestPerClonotypeFrame:
    def test_has_all_required_columns(self):
        rows = [{
            "clonotypeKey": "A", "kd": 10.0, "hillCoefficient": 1.0,
            "r2": 0.95, "affinityClass": "Good", "fitFailureReason": None,
            "kdOutOfRange": False,
        }]
        df = build_per_clonotype_frame(rows)
        assert set(df.columns) == {
            "clonotypeKey", "kd", "hillCoefficient", "r2",
            "affinityClass", "fitFailureReason", "kdOutOfRange",
        }

    def test_failed_row_has_null_kd(self):
        rows = [{
            "clonotypeKey": "A", "kd": None, "hillCoefficient": None,
            "r2": None, "affinityClass": "Failed", "fitFailureReason": "low_r2",
            "kdOutOfRange": None,
        }]
        df = build_per_clonotype_frame(rows)
        assert df["kd"][0] is None
        assert df["fitFailureReason"][0] == "low_r2"

    def test_empty_rows_preserves_schema(self):
        df = build_per_clonotype_frame([])
        assert df.height == 0
        assert "kd" in df.columns


class TestFlagKdOutOfRange:
    # R14b: outside [min, max] → kdOutOfRange = True. Closed interval (boundary = in-range).
    def test_in_range_false(self):
        df = build_per_clonotype_frame([{
            "clonotypeKey": "A", "kd": 1.0, "hillCoefficient": 1.0, "r2": 0.9,
            "affinityClass": "Good", "fitFailureReason": None, "kdOutOfRange": None,
        }])
        out = flag_kd_out_of_range(df, min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is False

    def test_above_max_flag_true(self):
        df = build_per_clonotype_frame([{
            "clonotypeKey": "A", "kd": 1000.0, "hillCoefficient": 1.0, "r2": 0.9,
            "affinityClass": "Good", "fitFailureReason": None, "kdOutOfRange": None,
        }])
        out = flag_kd_out_of_range(df, min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is True

    def test_below_min_flag_true(self):
        df = build_per_clonotype_frame([{
            "clonotypeKey": "A", "kd": 0.01, "hillCoefficient": 1.0, "r2": 0.9,
            "affinityClass": "Good", "fitFailureReason": None, "kdOutOfRange": None,
        }])
        out = flag_kd_out_of_range(df, min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is True

    @pytest.mark.parametrize("kd_value", [0.1, 100.0])
    def test_exactly_at_boundary_in_range(self, kd_value):
        df = build_per_clonotype_frame([{
            "clonotypeKey": "A", "kd": kd_value, "hillCoefficient": 1.0, "r2": 0.9,
            "affinityClass": "Good", "fitFailureReason": None, "kdOutOfRange": None,
        }])
        out = flag_kd_out_of_range(df, min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is False

    def test_failed_fit_kd_null_stays_null(self):
        df = build_per_clonotype_frame([{
            "clonotypeKey": "A", "kd": None, "hillCoefficient": None, "r2": None,
            "affinityClass": "Failed", "fitFailureReason": "low_r2", "kdOutOfRange": None,
        }])
        out = flag_kd_out_of_range(df, min_concentration=0.1, max_concentration=100.0)
        assert out["kdOutOfRange"][0] is None


class TestMeanBinFrame:
    def test_c0_excluded(self):
        signal = pl.DataFrame([
            {"clonotypeKey": "A", "concentrationStr": "0", "concentration": 0.0, "signal": 1.2},
            {"clonotypeKey": "A", "concentrationStr": "10", "concentration": 10.0, "signal": 2.5},
        ])
        out = build_mean_bin_frame(signal)
        assert out.height == 1
        assert out["concentrationStr"][0] == "10"
        assert out["meanBin"][0] == 2.5

    def test_preserves_canonical_string(self):
        signal = pl.DataFrame([
            {"clonotypeKey": "A", "concentrationStr": "1.000", "concentration": 1.0, "signal": 2.0},
        ])
        out = build_mean_bin_frame(signal)
        assert out["concentrationStr"][0] == "1.000"


class TestFittedMeanBinFrame:
    def test_empty_rows_preserves_schema(self):
        df = build_fitted_mean_bin_frame([])
        assert df.height == 0
        assert set(df.columns) == {"clonotypeKey", "concentrationStr", "concentration", "fittedMeanBin"}

    def test_populated_rows(self):
        rows = [{"clonotypeKey": "A", "concentrationStr": "10", "concentration": 10.0, "fittedMeanBin": 2.7}]
        df = build_fitted_mean_bin_frame(rows)
        assert df["fittedMeanBin"][0] == pytest.approx(2.7)
