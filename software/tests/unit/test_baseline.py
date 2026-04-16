"""Behavioral tests for baseline.py (R6 global baseline B)."""

from __future__ import annotations

import polars as pl
import pytest

from baseline import compute_global_baseline


def _c0(signals):
    return pl.DataFrame({
        "clonotypeKey": [f"C{i}" for i in range(len(signals))],
        "concentrationStr": ["0"] * len(signals),
        "signal": signals,
        "weight": [10.0] * len(signals),
    })


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
            schema={"clonotypeKey": pl.Utf8, "concentrationStr": pl.Utf8,
                    "signal": pl.Float64, "weight": pl.Float64},
        )
        assert compute_global_baseline(df) is None

    # Upstream floor filter drops some clonotypes from the c=0 set;
    # only surviving rows contribute to B.
    def test_only_passes_surviving_rows(self):
        df = _c0([1.2, 1.8])  # pretend two survived; floor-filtered third absent
        assert compute_global_baseline(df) == pytest.approx(1.5)
