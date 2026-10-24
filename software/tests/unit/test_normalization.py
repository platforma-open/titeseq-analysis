"""Behavioral tests for normalization.py (R7 mean bin, R7b frequency signal)."""

from __future__ import annotations

import polars as pl
import pytest

from normalization import MEAN_BIN, SIGNAL, compute_frequency_signal, compute_mean_bin, normalize


def _read(
    clonotype: str,
    reads: int,
    *,
    sample: str = "s",
    conc_str: str = "1",
    bin_label: int | None = None,
    sort_fraction: float | None = None,
) -> dict:
    """One reads row. The bin and sort_fraction columns are added only when given."""
    row: dict = {
        "clonotypeKey": clonotype,
        "sampleId": sample,
        "concentrationStr": conc_str,
        "concentration": float(conc_str),
    }
    if bin_label is not None:
        row["bin"] = bin_label
    row["reads"] = reads
    if sort_fraction is not None:
        row["sort_fraction"] = sort_fraction
    return row


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
        rows.append(_read("A", r, sample=f"s_b{bin_label}", bin_label=bin_label, sort_fraction=frac))
        if d - r > 0:
            rows.append(_read("F", d - r, sample=f"s_b{bin_label}", bin_label=bin_label, sort_fraction=frac))
    return pl.DataFrame(rows)


class TestComputeMeanBin:
    # Mean_bin_c = Σ(b·freq_cb) / Σ(freq_cb). Guards against using raw counts.
    @pytest.mark.parametrize(
        "reads_per_bin, depth_per_bin, expected_mean_bin",
        [
            # freqs = [0.01, 0.04, 0.025, 0.05] → 0.365 / 0.125 = 2.92
            ([10, 20, 5, 5], [1000, 500, 200, 100], pytest.approx(2.92, abs=1e-9)),
            # Single bin has reads → mean_bin == that bin label; zero-read bins add to neither sum
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
        rows = [_read("A", reads_clonotype, conc_str="10")]
        if reads_other > 0:
            rows.append(_read("B", reads_other, conc_str="10"))
        reads = pl.DataFrame(rows)
        out = compute_frequency_signal(reads).filter(pl.col("clonotypeKey") == "A")
        assert out.height == 1
        assert out[SIGNAL][0] == pytest.approx(expected)

    # total_reads_at_conc = 0 → signal is null (row survives with null, caller drops).
    def test_zero_total_reads_yields_null_signal(self):
        reads = pl.DataFrame([_read("A", 0, conc_str="10")])
        out = compute_frequency_signal(reads)
        assert out[SIGNAL][0] is None


class TestComputeMeanBinSortFraction:
    """FACS sort-fraction correction (Adams, Mora, Walczak, Kinney 2016 eq. A3).

    `sort_fraction_col=None` is the regression gate — must stay bit-exact with
    legacy behaviour. Non-None activates the weighted-mean correction.
    """

    # Corrected mean_bin = Σ(b·freq·frac) / Σ(freq·frac), freq = reads/depth per bin.
    @pytest.mark.parametrize(
        "reads_per_bin, depth_per_bin, fractions, expected",
        [
            # Equal fractions make the correction a constant factor: identical to the uncorrected 2.92.
            pytest.param([10, 20, 5, 5], [1000, 500, 200, 100], [0.25] * 4, 2.92, id="uniform-equals-uncorrected"),
            # freq=[0.1, 0.1], frac=[0.8, 0.2]: (1·0.08 + 2·0.02) / 0.10 = 1.2, pulled toward the heavy bin.
            pytest.param([10, 10], [100, 100], [0.8, 0.2], 1.2, id="skew-pulls-toward-heavy-bin"),
            # frac=[0.5, 0, 0.5]: bin 2's reads drop out of both sums, (0.05 + 0.15) / 0.10 = 2.0.
            pytest.param([10, 10, 10], [100, 100, 100], [0.5, 0.0, 0.5], 2.0, id="zero-fraction-bin-excluded"),
        ],
    )
    def test_corrected_mean_bin_pinned(self, reads_per_bin, depth_per_bin, fractions, expected):
        reads = _build_single_clonotype(reads_per_bin, depth_per_bin, sort_fraction_per_bin=fractions)
        out = compute_mean_bin(reads, sort_fraction_col="sort_fraction").filter(pl.col("clonotypeKey") == "A")
        assert out["mean_bin"][0] == pytest.approx(expected, abs=1e-9)

    # The correction must act inside each (clonotype, concentration) group.
    # A bug where the weighting leaks across concentrations would produce the
    # same `mean_bin` at both concs whenever the clonotype's reads are the same.
    # Use two concentrations with contrasting sort skews; they must produce
    # different mean_bin values.
    def test_per_concentration_independence(self):
        rows = []
        # conc 1: bin 1 dominates the sort (0.8, 0.2) → corrected mean_bin pulled toward 1.
        # conc 2: bin 2 dominates the sort (0.2, 0.8) → corrected mean_bin pulled toward 2.
        for conc_str, fractions in [("1", [0.8, 0.2]), ("2", [0.2, 0.8])]:
            for j, frac in enumerate(fractions):
                sample = f"s_c{conc_str}_b{j + 1}"
                for clonotype, n in [("A", 10), ("F", 90)]:
                    rows.append(
                        _read(clonotype, n, sample=sample, conc_str=conc_str, bin_label=j + 1, sort_fraction=frac)
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

    # Bin mode is compute_mean_bin with the column renamed; in no-bin mode sort_fraction (bin-space) is a no-op.
    def test_dispatch_by_mode(self):
        reads = _build_single_clonotype([10, 20, 5, 5], [1000, 500, 200, 100])
        keys = ["clonotypeKey", "concentrationStr"]
        direct = compute_mean_bin(reads, sort_fraction_col=None).rename({MEAN_BIN: SIGNAL}).sort(keys)
        assert normalize(reads, bin_mode=True, sort_fraction_col=None).sort(keys).equals(direct)
        no_bin = pl.DataFrame([_read("A", 100), _read("B", 900)])
        without = normalize(no_bin, bin_mode=False, sort_fraction_col=None).sort("clonotypeKey")
        assert normalize(no_bin, bin_mode=False, sort_fraction_col="sort_fraction").sort("clonotypeKey").equals(without)
