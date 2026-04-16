"""Behavioral tests for hook_effect.py (R9b non-monotonic signal detection)."""

from __future__ import annotations

import polars as pl
import pytest

from constants import DEFAULT_PARAMS, FitParams
from hook_effect import detect_hook_effect


def _build_fit_points(top2_signal, top1_signal, top2_reads, top1_reads):
    """One clonotype with two concentration points; 10 and 100 are the top-2 and top-1."""
    return pl.DataFrame(
        [
            {"clonotypeKey": "A", "concentrationStr": "10", "concentration": 10.0,
             "signal": top2_signal, "clonotype_reads_at_conc": top2_reads, "weight": float(top2_reads)},
            {"clonotypeKey": "A", "concentrationStr": "100", "concentration": 100.0,
             "signal": top1_signal, "clonotype_reads_at_conc": top1_reads, "weight": float(top1_reads)},
        ]
    )


class TestHookEffectBinMode:
    """threshold = 0.2; min_reads = 20. Strict > drop flags a hook."""

    @pytest.mark.parametrize(
        "top2_signal, top1_signal, top2_reads, top1_reads, expected",
        [
            (3.0, 2.5, 100, 100, True),    # drop 0.5 > 0.2 → flag
            (3.0, 2.85, 100, 100, False),  # drop 0.15 < 0.2 → no flag
            (3.0, 3.2, 100, 100, False),   # signal rose → no flag
            (3.0, 2.0, 100, 10, False),    # top1_reads < 20 → skip
            (3.0, 2.0, 10, 100, False),    # top2_reads < 20 → skip
            (3.0, 3.0, 100, 100, False),   # flat (zero drop) → no flag (strict >)
            (3.0, 2.81, 100, 100, False),  # drop 0.19 just under threshold → no flag
            (3.0, 2.79, 100, 100, True),   # drop 0.21 just over threshold → flag
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
            (0.05, 0.02, True),    # drop 0.03 > 0.02 → flag
            (0.05, 0.04, False),   # drop 0.01 < 0.02 → no flag
            (0.05, 0.031, False),   # drop 0.019 just under threshold → no flag
            (0.05, 0.029, True),    # drop 0.021 just over threshold → flag
        ],
    )
    def test_no_bin_mode_hook_detection(self, top2_signal, top1_signal, expected):
        df = _build_fit_points(top2_signal, top1_signal, 100, 100)
        result = detect_hook_effect(df, bin_mode=False, params=DEFAULT_PARAMS)
        row = result.filter(pl.col("clonotypeKey") == "A")
        assert bool(row["hook_flag"][0]) is expected
