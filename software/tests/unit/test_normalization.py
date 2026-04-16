"""Behavioral tests for normalization.py (R7 mean bin, R7b frequency signal)."""

from __future__ import annotations

import polars as pl
import pytest

from normalization import SIGNAL, compute_frequency_signal, compute_mean_bin


def _build_single_clonotype(reads_per_bin: list[int], depth_per_bin: list[int]) -> pl.DataFrame:
    """One clonotype at one concentration; multiple bins. Each bin's depth is set by adding
    filler reads under a dummy second clonotype so sum(reads at (bin, conc)) == depth."""
    rows = []
    for j, (r, d) in enumerate(zip(reads_per_bin, depth_per_bin)):
        bin_label = j + 1
        rows.append(
            {"clonotypeKey": "A", "sampleId": f"s_b{bin_label}",
             "concentrationStr": "1", "concentration": 1.0, "bin": bin_label, "reads": r}
        )
        filler = d - r
        if filler > 0:
            rows.append(
                {"clonotypeKey": "F", "sampleId": f"s_b{bin_label}",
                 "concentrationStr": "1", "concentration": 1.0, "bin": bin_label, "reads": filler}
            )
    return pl.DataFrame(rows)


class TestComputeMeanBin:
    # Mean_bin_c = Σ(b·freq_cb) / Σ(freq_cb). Guards against using raw counts.
    @pytest.mark.parametrize(
        "reads_per_bin, depth_per_bin, expected_mean_bin",
        [
            # freqs = [0.01, 0.04, 0.025, 0.05] → 0.365 / 0.125 = 2.92
            ([10, 20, 5, 5], [1000, 500, 200, 100], pytest.approx(2.92, abs=1e-9)),
            # Single bin has reads → mean_bin == that bin label
            ([0, 0, 7, 0], [1000, 500, 200, 100], 3.0),
            # Uniform frequency across 4 bins → centre of mass
            ([100, 50, 20, 10], [100, 50, 20, 10], 2.5),
            # Ceiling: reads only in the max bin label
            ([0, 0, 0, 0, 0, 0, 0, 5], [100, 100, 100, 100, 100, 100, 100, 100], 8.0),
        ],
    )
    def test_mean_bin_pinned(self, reads_per_bin, depth_per_bin, expected_mean_bin):
        reads = _build_single_clonotype(reads_per_bin, depth_per_bin)
        out = compute_mean_bin(reads).filter(pl.col("clonotypeKey") == "A")
        assert out.height == 1
        assert out["mean_bin"][0] == expected_mean_bin

    # Zero reads in one bin: that bin contributes 0 to both numerator and denominator.
    def test_zero_reads_bin_excluded_from_contribution(self):
        reads = _build_single_clonotype([0, 10, 0, 0], [100, 100, 100, 100])
        out = compute_mean_bin(reads).filter(pl.col("clonotypeKey") == "A")
        assert out["mean_bin"][0] == 2.0

    # Single-bin clonotype: mean_bin equals that bin label (pass-through).
    def test_single_bin_equals_that_bin_label(self):
        reads = pl.DataFrame(
            [
                {"clonotypeKey": "A", "sampleId": "s", "concentrationStr": "1",
                 "concentration": 1.0, "bin": 4, "reads": 7},
                {"clonotypeKey": "F", "sampleId": "s", "concentrationStr": "1",
                 "concentration": 1.0, "bin": 4, "reads": 93},
            ]
        )
        out = compute_mean_bin(reads).filter(pl.col("clonotypeKey") == "A")
        assert out["mean_bin"][0] == 4.0


class TestComputeFrequencySignal:
    # R7b: signal = reads_clonotype / total_reads_at_conc.
    @pytest.mark.parametrize(
        "reads_clonotype, reads_other, expected",
        [
            (100, 9_900, 0.01),
            (1, 9_999, 0.0001),
            (10_000, 0, 1.0),   # ceiling: every read is this clonotype
            (0, 10_000, 0.0),    # zero numerator
        ],
    )
    def test_no_bin_signal_pinned(self, reads_clonotype, reads_other, expected):
        rows = [{
            "clonotypeKey": "A", "sampleId": "s", "concentrationStr": "10",
            "concentration": 10.0, "reads": reads_clonotype,
        }]
        if reads_other > 0:
            rows.append({
                "clonotypeKey": "B", "sampleId": "s", "concentrationStr": "10",
                "concentration": 10.0, "reads": reads_other,
            })
        reads = pl.DataFrame(rows)
        out = compute_frequency_signal(reads).filter(pl.col("clonotypeKey") == "A")
        assert out.height == 1
        assert out[SIGNAL][0] == pytest.approx(expected)

    # total_reads_at_conc = 0 → signal is null (row survives with null, caller drops).
    def test_zero_total_reads_yields_null_signal(self):
        reads = pl.DataFrame(
            [{"clonotypeKey": "A", "sampleId": "s", "concentrationStr": "10",
              "concentration": 10.0, "reads": 0}]
        )
        out = compute_frequency_signal(reads)
        assert out[SIGNAL][0] is None
