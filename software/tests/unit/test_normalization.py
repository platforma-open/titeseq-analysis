"""Behavioral tests for normalization.py (R7 mean bin, R7b frequency signal)."""

from __future__ import annotations

import polars as pl
import pytest

from normalization import MEAN_BIN, SIGNAL, compute_frequency_signal, compute_mean_bin, normalize


def _build_single_clonotype(
    reads_per_bin: list[int],
    depth_per_bin: list[int],
    sort_fraction_per_bin: list[float] | None = None,
) -> pl.DataFrame:
    """One clonotype at one concentration; multiple bins. Each bin's depth is set by adding
    filler reads under a dummy second clonotype so sum(reads at (bin, conc)) == depth.

    When `sort_fraction_per_bin` is supplied, every row (clonotype "A" plus filler
    "F") at bin `j` carries the same sort_fraction value — sort fractions live on
    the sample, not the clonotype.
    """
    rows = []
    for j, (r, d) in enumerate(zip(reads_per_bin, depth_per_bin)):
        bin_label = j + 1
        frac = sort_fraction_per_bin[j] if sort_fraction_per_bin is not None else None
        row_a = {
            "clonotypeKey": "A",
            "sampleId": f"s_b{bin_label}",
            "concentrationStr": "1",
            "concentration": 1.0,
            "bin": bin_label,
            "reads": r,
        }
        if frac is not None:
            row_a["sort_fraction"] = frac
        rows.append(row_a)
        filler = d - r
        if filler > 0:
            row_f = {
                "clonotypeKey": "F",
                "sampleId": f"s_b{bin_label}",
                "concentrationStr": "1",
                "concentration": 1.0,
                "bin": bin_label,
                "reads": filler,
            }
            if frac is not None:
                row_f["sort_fraction"] = frac
            rows.append(row_f)
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
                {
                    "clonotypeKey": "A",
                    "sampleId": "s",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 4,
                    "reads": 7,
                },
                {
                    "clonotypeKey": "F",
                    "sampleId": "s",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 4,
                    "reads": 93,
                },
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
            (10_000, 0, 1.0),  # ceiling: every read is this clonotype
            (0, 10_000, 0.0),  # zero numerator
        ],
    )
    def test_no_bin_signal_pinned(self, reads_clonotype, reads_other, expected):
        rows = [
            {
                "clonotypeKey": "A",
                "sampleId": "s",
                "concentrationStr": "10",
                "concentration": 10.0,
                "reads": reads_clonotype,
            }
        ]
        if reads_other > 0:
            rows.append(
                {
                    "clonotypeKey": "B",
                    "sampleId": "s",
                    "concentrationStr": "10",
                    "concentration": 10.0,
                    "reads": reads_other,
                }
            )
        reads = pl.DataFrame(rows)
        out = compute_frequency_signal(reads).filter(pl.col("clonotypeKey") == "A")
        assert out.height == 1
        assert out[SIGNAL][0] == pytest.approx(expected)

    # total_reads_at_conc = 0 → signal is null (row survives with null, caller drops).
    def test_zero_total_reads_yields_null_signal(self):
        reads = pl.DataFrame(
            [{"clonotypeKey": "A", "sampleId": "s", "concentrationStr": "10", "concentration": 10.0, "reads": 0}]
        )
        out = compute_frequency_signal(reads)
        assert out[SIGNAL][0] is None


class TestComputeMeanBinSortFraction:
    """FACS sort-fraction correction (Adams, Mora, Walczak, Kinney 2016 eq. A3).

    `sort_fraction_col=None` is the regression gate — must stay bit-exact with
    legacy behaviour. Non-None activates the weighted-mean correction.
    """

    # The legacy default must be preserved. The pre-FACS test_mean_bin_pinned
    # cases above already anchor the uncorrected output numerically; this test
    # is the explicit guard that passing sort_fraction_col=None is identical
    # to the (deprecated) positional-argument-less call.
    def test_legacy_mode_matches_uncorrected(self):
        reads = _build_single_clonotype([10, 20, 5, 5], [1000, 500, 200, 100])
        legacy = compute_mean_bin(reads, sort_fraction_col=None).filter(pl.col("clonotypeKey") == "A")
        assert legacy["mean_bin"][0] == pytest.approx(2.92, abs=1e-9)

    # When every bin catches an equal cell fraction, the correction factor is
    # constant across bins → reduces bit-exactly to the legacy formula.
    def test_uniform_fractions_match_legacy(self):
        reads_per_bin = [10, 20, 5, 5]
        depth_per_bin = [1000, 500, 200, 100]
        n_bins = len(reads_per_bin)
        uniform = [1.0 / n_bins] * n_bins
        reads = _build_single_clonotype(reads_per_bin, depth_per_bin, sort_fraction_per_bin=uniform)
        corrected = compute_mean_bin(reads, sort_fraction_col="sort_fraction").filter(pl.col("clonotypeKey") == "A")
        legacy = compute_mean_bin(reads, sort_fraction_col=None).filter(pl.col("clonotypeKey") == "A")
        assert corrected["mean_bin"][0] == pytest.approx(legacy["mean_bin"][0], abs=1e-12)

    # Anchored 2-bin case: freq=[0.1, 0.1], fraction=[0.8, 0.2].
    # Uncorrected mean_bin = (1·0.1 + 2·0.1) / (0.1 + 0.1) = 1.5.
    # Corrected num = 1·(0.1·0.8) + 2·(0.1·0.2) = 0.12.
    # Corrected den = 0.1·0.8 + 0.1·0.2 = 0.10.
    # Corrected mean_bin = 1.2.
    def test_skewed_fractions_pinned_shifts_mean_toward_heavy_bin(self):
        reads = _build_single_clonotype(
            reads_per_bin=[10, 10],
            depth_per_bin=[100, 100],
            sort_fraction_per_bin=[0.8, 0.2],
        )
        out = compute_mean_bin(reads, sort_fraction_col="sort_fraction").filter(pl.col("clonotypeKey") == "A")
        assert out["mean_bin"][0] == pytest.approx(1.2)

    # Exclusion behaviour: a bin with sort_fraction=0 contributes to neither
    # numerator nor denominator even when it carries nonzero reads. This is the
    # formula acting correctly — multiplying freq by 0 zeros both num and den
    # contributions.
    def test_zero_sort_fraction_bin_excluded(self):
        # freq=[0.1, 0.1, 0.1] with fraction=[0.5, 0.0, 0.5].
        # Weighted num = 1·(0.1·0.5) + 2·(0.1·0.0) + 3·(0.1·0.5) = 0.05 + 0 + 0.15 = 0.20.
        # Weighted den = 0.1·0.5 + 0.1·0.0 + 0.1·0.5 = 0.10.
        # Corrected mean_bin = 2.0 — bin 2's reads silently removed, exactly as the spec intends.
        reads = _build_single_clonotype(
            reads_per_bin=[10, 10, 10],
            depth_per_bin=[100, 100, 100],
            sort_fraction_per_bin=[0.5, 0.0, 0.5],
        )
        out = compute_mean_bin(reads, sort_fraction_col="sort_fraction").filter(pl.col("clonotypeKey") == "A")
        assert out["mean_bin"][0] == pytest.approx(2.0)

    # The correction must act inside each (clonotype, concentration) group.
    # A bug where the weighting leaks across concentrations would produce the
    # same `mean_bin` at both concs whenever the clonotype's reads are the same.
    # Use two concentrations with contrasting sort skews; they must produce
    # different mean_bin values.
    def test_per_concentration_independence(self):
        rows = []
        # conc 1: bin 1 dominates the sort (0.8, 0.2) → corrected mean_bin pulled toward 1.
        # conc 2: bin 2 dominates the sort (0.2, 0.8) → corrected mean_bin pulled toward 2.
        for conc_str, conc_val, fractions in [("1", 1.0, [0.8, 0.2]), ("2", 2.0, [0.2, 0.8])]:
            for j, frac in enumerate(fractions):
                bin_label = j + 1
                rows.append(
                    {
                        "clonotypeKey": "A",
                        "sampleId": f"s_c{conc_str}_b{bin_label}",
                        "concentrationStr": conc_str,
                        "concentration": conc_val,
                        "bin": bin_label,
                        "reads": 10,
                        "sort_fraction": frac,
                    }
                )
                rows.append(
                    {
                        "clonotypeKey": "F",
                        "sampleId": f"s_c{conc_str}_b{bin_label}",
                        "concentrationStr": conc_str,
                        "concentration": conc_val,
                        "bin": bin_label,
                        "reads": 90,
                        "sort_fraction": frac,
                    }
                )
        reads = pl.DataFrame(rows)
        out = (
            compute_mean_bin(reads, sort_fraction_col="sort_fraction")
            .filter(pl.col("clonotypeKey") == "A")
            .sort("concentrationStr")
        )
        # Corrected mean_bin at conc 1: (1·0.1·0.8 + 2·0.1·0.2) / (0.1·0.8 + 0.1·0.2) = 0.12/0.10 = 1.2.
        # Corrected mean_bin at conc 2: (1·0.1·0.2 + 2·0.1·0.8) / (0.1·0.2 + 0.1·0.8) = 0.18/0.10 = 1.8.
        vals = {row["concentrationStr"]: row["mean_bin"] for row in out.to_dicts()}
        assert vals["1"] == pytest.approx(1.2)
        assert vals["2"] == pytest.approx(1.8)


class TestNormalizeDispatch:
    """normalize() threads sort_fraction_col to compute_mean_bin without mutating legacy paths."""

    def test_legacy_bin_mode_bit_exact(self):
        reads = _build_single_clonotype([10, 20, 5, 5], [1000, 500, 200, 100])
        legacy = normalize(reads, bin_mode=True, sort_fraction_col=None)
        direct = compute_mean_bin(reads, sort_fraction_col=None).rename({MEAN_BIN: SIGNAL})
        assert legacy.equals(direct)

    def test_frequency_mode_ignores_sort_fraction(self):
        # sort_fraction is bin-space; in no-bin mode the kwarg must be a no-op.
        # compute_frequency_signal never sees the column, so passing the kwarg
        # simply dispatches to the frequency path unchanged.
        reads = pl.DataFrame(
            [
                {"clonotypeKey": "A", "sampleId": "s", "concentrationStr": "1", "concentration": 1.0, "reads": 100},
                {"clonotypeKey": "B", "sampleId": "s", "concentrationStr": "1", "concentration": 1.0, "reads": 900},
            ]
        )
        out_without = normalize(reads, bin_mode=False, sort_fraction_col=None)
        out_with = normalize(reads, bin_mode=False, sort_fraction_col="sort_fraction")
        assert out_without.equals(out_with)
