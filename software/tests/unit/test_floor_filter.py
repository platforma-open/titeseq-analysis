"""Behavioral tests for floor_filter.py (R8 weights, R9 floor + insufficient-* classification)."""

from __future__ import annotations

import polars as pl
import pytest

from constants import DEFAULT_PARAMS, FitParams
from floor_filter import (
    WEIGHT,
    apply_floor_and_weights,
    classify_insufficient,
    split_c0,
)


def _sig_frame(rows):
    """Build a signal-frame shape: clonotypeKey, concentrationStr, concentration, signal, clonotype_reads_at_conc."""
    return pl.DataFrame(rows)


class TestApplyFloorAndWeights:
    # R8: weight w_j = clonotype_reads_at_conc at concentration j.
    def test_weight_equals_reads_sum(self):
        df = _sig_frame([
            {"clonotypeKey": "A", "concentrationStr": "1", "concentration": 1.0,
             "signal": 2.5, "clonotype_reads_at_conc": 10},
        ])
        result = apply_floor_and_weights(df, DEFAULT_PARAMS)
        assert result[WEIGHT][0] == 10.0

    # R9 floor: inclusive at the threshold (>= min_reads_per_concentration).
    # Below → dropped; at → kept; above → kept.
    @pytest.mark.parametrize(
        "reads_at_conc, kept",
        [(2, False), (3, True), (10, True)],
    )
    def test_floor_threshold_inclusive(self, reads_at_conc, kept):
        params = FitParams(min_reads_per_concentration=3)
        df = _sig_frame([
            {"clonotypeKey": "A", "concentrationStr": "1", "concentration": 1.0,
             "signal": 2.5, "clonotype_reads_at_conc": reads_at_conc},
        ])
        result = apply_floor_and_weights(df, params)
        assert (result.height == 1) is kept


class TestClassifyInsufficient:
    # Spec R9: zero surviving points → insufficient_reads.
    def test_zero_points_marks_insufficient_reads(self):
        filtered = pl.DataFrame(
            {"clonotypeKey": [], "concentrationStr": [], "concentration": [],
             "signal": [], "clonotype_reads_at_conc": [], WEIGHT: []},
            schema={"clonotypeKey": pl.Utf8, "concentrationStr": pl.Utf8,
                    "concentration": pl.Float64, "signal": pl.Float64,
                    "clonotype_reads_at_conc": pl.Int64, WEIGHT: pl.Float64},
        )
        result = classify_insufficient(filtered, ["A"], DEFAULT_PARAMS)
        row = result.filter(pl.col("clonotypeKey") == "A")
        assert row["insufficient_reason"][0] == "insufficient_reads"

    # R9: fewer than min_concentration_points → insufficient_points.
    def test_under_min_points_marks_insufficient_points(self):
        params = FitParams(min_concentration_points=5)
        rows = [
            {"clonotypeKey": "A", "concentrationStr": str(c), "concentration": float(c),
             "signal": 2.0, "clonotype_reads_at_conc": 10, WEIGHT: 10.0}
            for c in [1, 10, 100]  # only 3 non-zero points
        ]
        filtered = pl.DataFrame(rows)
        result = classify_insufficient(filtered, ["A"], params)
        assert result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"][0] == "insufficient_points"

    # Boundary: exactly min_concentration_points → not marked.
    def test_exactly_min_points_not_marked(self):
        params = FitParams(min_concentration_points=5)
        rows = [
            {"clonotypeKey": "A", "concentrationStr": str(c), "concentration": float(c),
             "signal": 2.0, "clonotype_reads_at_conc": 10, WEIGHT: 10.0}
            for c in [1, 3, 10, 30, 100]  # exactly 5
        ]
        filtered = pl.DataFrame(rows)
        result = classify_insufficient(filtered, ["A"], params)
        reason_cell = result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"]
        assert reason_cell.is_null().item() is True

    # c=0 rows are NOT counted toward the minimum points threshold (only fit points).
    def test_c0_excluded_from_point_count(self):
        params = FitParams(min_concentration_points=5)
        rows = [{"clonotypeKey": "A", "concentrationStr": "0", "concentration": 0.0,
                 "signal": 1.0, "clonotype_reads_at_conc": 10, WEIGHT: 10.0}]
        rows += [
            {"clonotypeKey": "A", "concentrationStr": str(c), "concentration": float(c),
             "signal": 2.0, "clonotype_reads_at_conc": 10, WEIGHT: 10.0}
            for c in [1, 3, 10, 30]  # only 4 non-zero → still insufficient
        ]
        filtered = pl.DataFrame(rows)
        result = classify_insufficient(filtered, ["A"], params)
        assert result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"][0] == "insufficient_points"


class TestSplitC0:
    def test_partitions_c0_from_fit_points(self):
        df = pl.DataFrame([
            {"clonotypeKey": "A", "concentrationStr": "0", "concentration": 0.0,
             "signal": 1.2, "clonotype_reads_at_conc": 100, WEIGHT: 100.0},
            {"clonotypeKey": "A", "concentrationStr": "10", "concentration": 10.0,
             "signal": 2.5, "clonotype_reads_at_conc": 100, WEIGHT: 100.0},
        ])
        c0, non_c0 = split_c0(df)
        assert c0.height == 1
        assert c0["concentrationStr"][0] == "0"
        assert non_c0.height == 1
