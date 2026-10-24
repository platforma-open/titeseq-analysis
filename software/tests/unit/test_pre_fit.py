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


def _c0(signals: list[float]) -> pl.DataFrame:
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

    # Upstream floor filter drops some clonotypes from the c=0 set;
    # only surviving rows contribute to B.
    def test_only_passes_surviving_rows(self):
        df = _c0([1.2, 1.8])  # pretend two survived; floor-filtered third absent
        assert compute_global_baseline(df) == pytest.approx(1.5)


def _point(conc_str: str, signal: float = 2.0, reads: int = 10, *, weighted: bool = True) -> dict:
    """One signal-frame row for clonotype A. `weighted` adds the weight column, equal to the reads."""
    row = {
        "clonotypeKey": "A",
        "concentrationStr": conc_str,
        "concentration": float(conc_str),
        "signal": signal,
        "clonotype_reads_at_conc": reads,
    }
    if weighted:
        row[WEIGHT] = float(reads)
    return row


class TestApplyFloorAndWeights:
    # R8: weight w_j = clonotype_reads_at_conc at concentration j.
    def test_weight_equals_reads_sum(self):
        df = pl.DataFrame([_point("1", 2.5, 10, weighted=False)])
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
        df = pl.DataFrame([_point("1", 2.5, reads_at_conc, weighted=False)])
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
        filtered = pl.DataFrame([_point(str(c)) for c in [1, 10, 100]])  # only 3 non-zero points
        result = classify_insufficient(filtered, ["A"], params)
        assert result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"][0] == "insufficient_points"

    # Boundary: exactly min_concentration_points → not marked.
    def test_exactly_min_points_not_marked(self):
        params = FitParams(min_concentration_points=5)
        filtered = pl.DataFrame([_point(str(c)) for c in [1, 3, 10, 30, 100]])  # exactly 5
        result = classify_insufficient(filtered, ["A"], params)
        reason_cell = result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"]
        assert reason_cell.is_null().item() is True

    # c=0 rows are NOT counted toward the minimum points threshold (only fit points).
    def test_c0_excluded_from_point_count(self):
        params = FitParams(min_concentration_points=5)
        # only 4 non-zero → still insufficient
        filtered = pl.DataFrame([_point("0", 1.0)] + [_point(str(c)) for c in [1, 3, 10, 30]])
        result = classify_insufficient(filtered, ["A"], params)
        assert result.filter(pl.col("clonotypeKey") == "A")["insufficient_reason"][0] == "insufficient_points"


class TestSplitC0:
    def test_partitions_c0_from_fit_points(self):
        df = pl.DataFrame([_point("0", 1.2, 100), _point("10", 2.5, 100)])
        c0, non_c0 = split_c0(df)
        assert c0.height == 1
        assert c0["concentrationStr"][0] == "0"
        assert non_c0.height == 1


def _fit_points(signals: tuple[float, float, float], reads: tuple[int, int, int] = (100, 100, 100)) -> pl.DataFrame:
    """One clonotype at concentrations 1 / 10 / 100: (top-3, top-2, top-1) signals and reads."""
    return pl.DataFrame([_point(c, signal, n) for c, signal, n in zip(("1", "10", "100"), signals, reads)])


# R9b: flag a hook when top2 - top1 > δ AND top3 - top1 > δ/2 (strict), with δ = 0.2 in bin mode and
# 0.02 in no-bin mode. Points with fewer than 20 reads at top-1 or top-2 are not evaluated; top-3 is not gated.
# Cases with top3 == top2 exercise the first clause alone.
@pytest.mark.parametrize(
    "bin_mode, signals, reads, expected",
    [
        pytest.param(True, (3.0, 3.0, 2.5), (100, 100, 100), True, id="bin-drop-0.5-flags"),
        pytest.param(True, (3.0, 3.0, 2.85), (100, 100, 100), False, id="bin-drop-0.15-under"),
        pytest.param(True, (3.0, 3.0, 3.2), (100, 100, 100), False, id="bin-signal-rose"),
        pytest.param(True, (3.0, 3.0, 3.0), (100, 100, 100), False, id="bin-flat-strict"),
        pytest.param(True, (3.0, 3.0, 2.81), (100, 100, 100), False, id="bin-drop-0.19-just-under"),
        pytest.param(True, (3.0, 3.0, 2.79), (100, 100, 100), True, id="bin-drop-0.21-just-over"),
        pytest.param(True, (2.52, 3.0, 2.5), (100, 100, 100), False, id="bin-top3-clause-fails"),
        pytest.param(True, (2.8, 3.0, 2.5), (100, 100, 100), True, id="bin-both-clauses"),
        pytest.param(True, (3.0, 2.95, 2.8), (100, 100, 100), False, id="bin-first-clause-fails"),
        pytest.param(True, (2.59, 3.0, 2.5), (100, 100, 100), False, id="bin-half-threshold-just-under"),
        pytest.param(True, (2.601, 3.0, 2.5), (100, 100, 100), True, id="bin-half-threshold-just-over"),
        pytest.param(True, (3.0, 3.0, 2.0), (100, 100, 10), False, id="bin-top1-below-read-floor"),
        pytest.param(True, (3.0, 3.0, 2.0), (10, 10, 100), False, id="bin-top2-below-read-floor"),
        pytest.param(True, (2.8, 3.0, 2.5), (10, 100, 100), True, id="bin-top3-reads-not-gated"),
        pytest.param(True, (2.8, 3.0, 2.5), (100, 10, 100), False, id="bin-top2-gated"),
        pytest.param(True, (2.8, 3.0, 2.5), (100, 100, 10), False, id="bin-top1-gated"),
        pytest.param(False, (0.05, 0.05, 0.02), (100, 100, 100), True, id="nobin-drop-0.03-flags"),
        pytest.param(False, (0.05, 0.05, 0.04), (100, 100, 100), False, id="nobin-drop-0.01-under"),
        pytest.param(False, (0.05, 0.05, 0.031), (100, 100, 100), False, id="nobin-drop-0.019-just-under"),
        pytest.param(False, (0.05, 0.05, 0.029), (100, 100, 100), True, id="nobin-drop-0.021-just-over"),
        pytest.param(False, (0.04, 0.05, 0.02), (100, 100, 100), True, id="nobin-both-clauses"),
        pytest.param(False, (0.025, 0.05, 0.02), (100, 100, 100), False, id="nobin-top3-clause-fails"),
        pytest.param(False, (0.029, 0.05, 0.02), (100, 100, 100), False, id="nobin-half-threshold-just-under"),
    ],
)
def test_hook_effect_detection(bin_mode, signals, reads, expected):
    result = detect_hook_effect(_fit_points(signals, reads), bin_mode=bin_mode, params=DEFAULT_PARAMS)
    assert bool(result.filter(pl.col("clonotypeKey") == "A")["hook_flag"][0]) is expected
