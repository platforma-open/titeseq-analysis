"""Behavioral tests for pre_fit.py (R6 baseline, R8 weights, R9 floor, R9b hook effect)."""

from __future__ import annotations

import polars as pl
import pytest

from constants import DEFAULT_PARAMS, FitParams
from pre_fit import (
    WEIGHT,
    apply_floor_and_weights,
    classify_insufficient,
    compute_global_baseline,
    detect_hook_effect,
    split_c0,
)


def _c0(signals):
    return pl.DataFrame(
        {
            "clonotypeKey": [f"C{i}" for i in range(len(signals))],
            "concentrationStr": ["0"] * len(signals),
            "signal": signals,
            "weight": [10.0] * len(signals),
        }
    )


class TestGlobalBaseline:
    # Arithmetic mean of per-clonotype mean_bin at c=0; spec explicit.
    def test_arithmetic_mean_of_c0_signals(self):
        df = _c0([1.2, 1.4, 1.6])
        assert compute_global_baseline(df) == pytest.approx(1.4)

    # A c=0 signal of exactly 0.0 is a valid observation — counts in the mean.
    def test_zero_signal_counts_in_mean(self):
        df = _c0([0.0, 2.0])
        assert compute_global_baseline(df) == pytest.approx(1.0)

    # No c=0 data at all → None (downstream uses 4-param fit).
    def test_empty_returns_none(self):
        df = pl.DataFrame(
            {"clonotypeKey": [], "concentrationStr": [], "signal": [], "weight": []},
            schema={"clonotypeKey": pl.Utf8, "concentrationStr": pl.Utf8, "signal": pl.Float64, "weight": pl.Float64},
        )
        assert compute_global_baseline(df) is None

    # Upstream floor filter drops some clonotypes from the c=0 set;
    # only surviving rows contribute to B.
    def test_only_passes_surviving_rows(self):
        df = _c0([1.2, 1.8])  # pretend two survived; floor-filtered third absent
        assert compute_global_baseline(df) == pytest.approx(1.5)


def _sig_frame(rows):
    """Signal-frame shape: clonotypeKey, concentrationStr, concentration, signal, clonotype_reads_at_conc."""
    return pl.DataFrame(rows)


class TestApplyFloorAndWeights:
    # R8: weight w_j = clonotype_reads_at_conc at concentration j.
    def test_weight_equals_reads_sum(self):
        df = _sig_frame(
            [
                {
                    "clonotypeKey": "A",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "signal": 2.5,
                    "clonotype_reads_at_conc": 10,
                },
            ]
        )
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
        df = _sig_frame(
            [
                {
                    "clonotypeKey": "A",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "signal": 2.5,
                    "clonotype_reads_at_conc": reads_at_conc,
                },
            ]
        )
        result = apply_floor_and_weights(df, params)
        assert (result.height == 1) is kept


class TestClassifyInsufficient:
    # Spec R9: zero surviving points → insufficient_reads.
    def test_zero_points_marks_insufficient_reads(self):
        filtered = pl.DataFrame(
            {
                "clonotypeKey": [],
                "concentrationStr": [],
                "concentration": [],
                "signal": [],
                "clonotype_reads_at_conc": [],
                WEIGHT: [],
            },
            schema={
                "clonotypeKey": pl.Utf8,
                "concentrationStr": pl.Utf8,
                "concentration": pl.Float64,
                "signal": pl.Float64,
                "clonotype_reads_at_conc": pl.Int64,
                WEIGHT: pl.Float64,
            },
        )
        result = classify_insufficient(filtered, ["A"], DEFAULT_PARAMS)
        row = result.filter(pl.col("clonotypeKey") == "A")
        assert row["insufficient_reason"][0] == "insufficient_reads"

    # R9: fewer than min_concentration_points → insufficient_points.
    def test_under_min_points_marks_insufficient_points(self):
        params = FitParams(min_concentration_points=5)
        rows = [
            {
                "clonotypeKey": "A",
                "concentrationStr": str(c),
                "concentration": float(c),
                "signal": 2.0,
                "clonotype_reads_at_conc": 10,
                WEIGHT: 10.0,
            }
            for c in [1, 10, 100]  # only 3 non-zero points
        ]
        filtered = pl.DataFrame(rows)
        result = classify_insufficient(filtered, ["A"], params)
        assert result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"][0] == "insufficient_points"

    # Boundary: exactly min_concentration_points → not marked.
    def test_exactly_min_points_not_marked(self):
        params = FitParams(min_concentration_points=5)
        rows = [
            {
                "clonotypeKey": "A",
                "concentrationStr": str(c),
                "concentration": float(c),
                "signal": 2.0,
                "clonotype_reads_at_conc": 10,
                WEIGHT: 10.0,
            }
            for c in [1, 3, 10, 30, 100]  # exactly 5
        ]
        filtered = pl.DataFrame(rows)
        result = classify_insufficient(filtered, ["A"], params)
        reason_cell = result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"]
        assert reason_cell.is_null().item() is True

    # c=0 rows are NOT counted toward the minimum points threshold (only fit points).
    def test_c0_excluded_from_point_count(self):
        params = FitParams(min_concentration_points=5)
        rows = [
            {
                "clonotypeKey": "A",
                "concentrationStr": "0",
                "concentration": 0.0,
                "signal": 1.0,
                "clonotype_reads_at_conc": 10,
                WEIGHT: 10.0,
            }
        ]
        rows += [
            {
                "clonotypeKey": "A",
                "concentrationStr": str(c),
                "concentration": float(c),
                "signal": 2.0,
                "clonotype_reads_at_conc": 10,
                WEIGHT: 10.0,
            }
            for c in [1, 3, 10, 30]  # only 4 non-zero → still insufficient
        ]
        filtered = pl.DataFrame(rows)
        result = classify_insufficient(filtered, ["A"], params)
        assert result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"][0] == "insufficient_points"


class TestSplitC0:
    def test_partitions_c0_from_fit_points(self):
        df = pl.DataFrame(
            [
                {
                    "clonotypeKey": "A",
                    "concentrationStr": "0",
                    "concentration": 0.0,
                    "signal": 1.2,
                    "clonotype_reads_at_conc": 100,
                    WEIGHT: 100.0,
                },
                {
                    "clonotypeKey": "A",
                    "concentrationStr": "10",
                    "concentration": 10.0,
                    "signal": 2.5,
                    "clonotype_reads_at_conc": 100,
                    WEIGHT: 100.0,
                },
            ]
        )
        c0, non_c0 = split_c0(df)
        assert c0.height == 1
        assert c0["concentrationStr"][0] == "0"
        assert non_c0.height == 1


def _build_fit_points(top2_signal, top1_signal, top2_reads, top1_reads):
    """One clonotype with two concentration points; 10 and 100 are the top-2 and top-1."""
    return pl.DataFrame(
        [
            {
                "clonotypeKey": "A",
                "concentrationStr": "10",
                "concentration": 10.0,
                "signal": top2_signal,
                "clonotype_reads_at_conc": top2_reads,
                "weight": float(top2_reads),
            },
            {
                "clonotypeKey": "A",
                "concentrationStr": "100",
                "concentration": 100.0,
                "signal": top1_signal,
                "clonotype_reads_at_conc": top1_reads,
                "weight": float(top1_reads),
            },
        ]
    )


class TestHookEffectBinMode:
    """threshold = 0.2; min_reads = 20. Strict > drop flags a hook."""

    @pytest.mark.parametrize(
        "top2_signal, top1_signal, top2_reads, top1_reads, expected",
        [
            (3.0, 2.5, 100, 100, True),  # drop 0.5 > 0.2 → flag
            (3.0, 2.85, 100, 100, False),  # drop 0.15 < 0.2 → no flag
            (3.0, 3.2, 100, 100, False),  # signal rose → no flag
            (3.0, 2.0, 100, 10, False),  # top1_reads < 20 → skip
            (3.0, 2.0, 10, 100, False),  # top2_reads < 20 → skip
            (3.0, 3.0, 100, 100, False),  # flat (zero drop) → no flag (strict >)
            (3.0, 2.81, 100, 100, False),  # drop 0.19 just under threshold → no flag
            (3.0, 2.79, 100, 100, True),  # drop 0.21 just over threshold → flag
        ],
    )
    def test_bin_mode_hook_detection(self, top2_signal, top1_signal, top2_reads, top1_reads, expected):
        df = _build_fit_points(top2_signal, top1_signal, top2_reads, top1_reads)
        result = detect_hook_effect(df, bin_mode=True, params=DEFAULT_PARAMS)
        row = result.filter(pl.col("clonotypeKey") == "A")
        assert bool(row["hook_flag"][0]) is expected


class TestHookEffectNoBinMode:
    """threshold = 0.02 (default no-bin)."""

    @pytest.mark.parametrize(
        "top2_signal, top1_signal, expected",
        [
            (0.05, 0.02, True),  # drop 0.03 > 0.02 → flag
            (0.05, 0.04, False),  # drop 0.01 < 0.02 → no flag
            (0.05, 0.031, False),  # drop 0.019 just under threshold → no flag
            (0.05, 0.029, True),  # drop 0.021 just over threshold → flag
        ],
    )
    def test_no_bin_mode_hook_detection(self, top2_signal, top1_signal, expected):
        df = _build_fit_points(top2_signal, top1_signal, 100, 100)
        result = detect_hook_effect(df, bin_mode=False, params=DEFAULT_PARAMS)
        row = result.filter(pl.col("clonotypeKey") == "A")
        assert bool(row["hook_flag"][0]) is expected


def _build_fit_points_3(
    top3_signal,
    top2_signal,
    top1_signal,
    top3_reads=100,
    top2_reads=100,
    top1_reads=100,
):
    """One clonotype with three concentration points at 1 / 10 / 100 (top-3 / top-2 / top-1)."""
    return pl.DataFrame(
        [
            {
                "clonotypeKey": "A",
                "concentrationStr": "1",
                "concentration": 1.0,
                "signal": top3_signal,
                "clonotype_reads_at_conc": top3_reads,
                "weight": float(top3_reads),
            },
            {
                "clonotypeKey": "A",
                "concentrationStr": "10",
                "concentration": 10.0,
                "signal": top2_signal,
                "clonotype_reads_at_conc": top2_reads,
                "weight": float(top2_reads),
            },
            {
                "clonotypeKey": "A",
                "concentrationStr": "100",
                "concentration": 100.0,
                "signal": top1_signal,
                "clonotype_reads_at_conc": top1_reads,
                "weight": float(top1_reads),
            },
        ]
    )


class TestHookEffectBinModeTop3:
    """R9b (bin mode, δ=0.2): flag only when BOTH top2-top1 > δ AND top3-top1 > δ/2."""

    # top2-top1 > δ alone is not enough — the top-3 clause guards against
    # genuine dose-response curves that dip only at the very top concentration.
    @pytest.mark.parametrize(
        "top3, top2, top1, expected, case",
        [
            # top2-top1 = 0.5 > 0.2 (first clause passes), top3-top1 = 0.02 ≤ 0.1 (second fails)
            (2.52, 3.0, 2.5, False, "first_cond_met_second_not_met"),
            # top2-top1 = 0.5 > 0.2, top3-top1 = 0.3 > 0.1 → flag
            (2.8, 3.0, 2.5, True, "both_cond_met"),
            # top2-top1 = 0.15 < 0.2 → no flag regardless of top-3
            (3.0, 2.95, 2.8, False, "first_cond_not_met"),
            # top3-top1 = exactly δ/2 = 0.1 → strict > means NOT flagged
            (2.6, 3.0, 2.5, False, "half_threshold_boundary_equal"),
            # top3-top1 = 0.101 just over δ/2 → flag
            (2.601, 3.0, 2.5, True, "half_threshold_boundary_just_over"),
        ],
        ids=lambda x: x if isinstance(x, str) else None,
    )
    def test_bin_mode_top3_clause(self, top3, top2, top1, expected, case):
        df = _build_fit_points_3(top3, top2, top1)
        result = detect_hook_effect(df, bin_mode=True, params=DEFAULT_PARAMS)
        assert bool(result.filter(pl.col("clonotypeKey") == "A")["hook_flag"][0]) is expected, case

    # Spec R9b: min-reads floor applies to the top-1, top-2, AND top-3 points.
    # If any of the three is below floor the check is skipped entirely.
    @pytest.mark.parametrize(
        "top3_reads, top2_reads, top1_reads, expected",
        [
            (100, 100, 100, True),  # all ≥ 20 → evaluate; drop qualifies under both clauses
            (10, 100, 100, False),  # top-3 below floor → skip
            (100, 10, 100, False),  # top-2 below floor → skip
            (100, 100, 10, False),  # top-1 below floor → skip
        ],
    )
    def test_min_reads_gate_extends_to_top3(self, top3_reads, top2_reads, top1_reads, expected):
        # Signals chosen so both spec conditions are satisfied when the gate allows.
        df = _build_fit_points_3(
            2.8,
            3.0,
            2.5,
            top3_reads=top3_reads,
            top2_reads=top2_reads,
            top1_reads=top1_reads,
        )
        result = detect_hook_effect(df, bin_mode=True, params=DEFAULT_PARAMS)
        assert bool(result.filter(pl.col("clonotypeKey") == "A")["hook_flag"][0]) is expected


class TestHookEffectNoBinModeTop3:
    """R9b (no-bin mode, δ=0.02): same two-clause logic at a smaller scale."""

    @pytest.mark.parametrize(
        "top3, top2, top1, expected, case",
        [
            # top2-top1 = 0.03 > 0.02 AND top3-top1 = 0.02 > 0.01 → flag
            (0.04, 0.05, 0.02, True, "both_cond_met"),
            # top2-top1 = 0.03 > 0.02 but top3-top1 = 0.005 < 0.01 → NOT flagged
            (0.025, 0.05, 0.02, False, "first_cond_met_second_not_met"),
            # top3-top1 exactly δ/2 = 0.01 → strict > → NOT flagged
            (0.03, 0.05, 0.02, False, "half_threshold_boundary_equal"),
        ],
        ids=lambda x: x if isinstance(x, str) else None,
    )
    def test_no_bin_mode_top3_clause(self, top3, top2, top1, expected, case):
        df = _build_fit_points_3(top3, top2, top1)
        result = detect_hook_effect(df, bin_mode=False, params=DEFAULT_PARAMS)
        assert bool(result.filter(pl.col("clonotypeKey") == "A")["hook_flag"][0]) is expected, case
