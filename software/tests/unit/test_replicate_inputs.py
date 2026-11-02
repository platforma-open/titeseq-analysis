"""The replicate column: what makes it usable, and what it changes about the grouping."""

from __future__ import annotations

import polars as pl
import pytest

from io_layer import InputValidationError, validate_replicate_column, validate_sort_fraction
from normalization import compute_frequency_signal, compute_mean_bin
from output_build import build_mean_bin_frame

REPS = ["replicate"]


def _samples(replicates: list[str | None]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "sampleId": [f"s{i}" for i in range(len(replicates))],
            "replicate": replicates,
        }
    )


class TestValidateReplicateColumn:
    def test_one_label_a_sample_is_usable(self):
        assert validate_replicate_column(_samples(["A", "A", "B"]), "replicate") == []

    def test_a_missing_column_warns(self):
        warnings = validate_replicate_column(_samples(["A"]).drop("replicate"), "replicate")
        assert len(warnings) == 1
        assert "one replicate" in warnings[0]

    def test_an_empty_value_warns(self):
        assert validate_replicate_column(_samples(["A", None]), "replicate") != []
        assert validate_replicate_column(_samples(["A", ""]), "replicate") != []

    def test_two_labels_on_one_sample_warn(self):
        frame = pl.DataFrame({"sampleId": ["s0", "s0"], "replicate": ["A", "B"]})
        warnings = validate_replicate_column(frame, "replicate")
        assert "one value per sample" in warnings[0]

    def test_an_integer_label_is_usable(self):
        frame = pl.DataFrame({"sampleId": ["s0", "s1"], "replicate": [1, 2]})
        assert validate_replicate_column(frame, "replicate") == []


def _sorts(replicates: list[str]) -> pl.DataFrame:
    """One concentration, four gates, each replicate's fractions summing to one on their own."""
    rows = []
    for replicate in replicates:
        for gate in range(4):
            rows.append(
                {
                    "sampleId": f"{replicate}_{gate}",
                    "concentrationStr": "1e-08",
                    "replicate": replicate,
                    "sort_fraction": 0.25,
                }
            )
    return pl.DataFrame(rows)


class TestSortFractionGrouping:
    def test_two_replicates_sum_to_two_when_pooled(self):
        with pytest.raises(InputValidationError, match="per concentration"):
            validate_sort_fraction(_sorts(["A", "B"]), "sort_fraction")

    def test_two_replicates_each_sum_to_one_within_the_replicate(self):
        validate_sort_fraction(_sorts(["A", "B"]), "sort_fraction", [*REPS, "concentrationStr"])

    def test_a_replicate_whose_fractions_do_not_sum_to_one_is_still_refused(self):
        frame = _sorts(["A", "B"]).with_columns(
            pl.when(pl.col("replicate") == "B").then(0.2).otherwise(pl.col("sort_fraction")).alias("sort_fraction")
        )
        with pytest.raises(InputValidationError, match="replicate"):
            validate_sort_fraction(frame, "sort_fraction", [*REPS, "concentrationStr"])


def _reads(clonotype: str, replicate: str, reads: list[int]) -> list[dict]:
    return [
        {
            "clonotypeKey": clonotype,
            "replicate": replicate,
            "concentrationStr": "1e-08",
            "concentration": 1e-8,
            "bin": gate + 1,
            "reads": count,
        }
        for gate, count in enumerate(reads)
    ]


class TestDepthGrouping:
    def test_a_replicate_sorted_at_another_depth_moves_a_pooled_mean_bin(self):
        """Replicate B's gate 2 holds a thousand times A's reads, so a pooled gate depth divides
        A's gate-2 reads by B's library and all but erases them."""
        frame = pl.DataFrame(_reads("v1", "A", [90, 10]) + _reads("v2", "B", [100, 10000]))
        pooled = compute_mean_bin(frame).filter(pl.col("clonotypeKey") == "v1")["mean_bin"].item()
        apart = compute_mean_bin(frame, replicate_keys=REPS).filter(pl.col("clonotypeKey") == "v1")["mean_bin"].item()
        assert pooled == pytest.approx(1.0021, abs=1e-4)
        assert apart == pytest.approx(1.5)

    def test_the_library_depth_of_a_frequency_run_is_taken_within_the_replicate(self):
        frame = pl.DataFrame(
            [
                {
                    "clonotypeKey": "v1",
                    "replicate": "A",
                    "concentrationStr": "1e-08",
                    "concentration": 1e-8,
                    "reads": 10,
                },
                {
                    "clonotypeKey": "v2",
                    "replicate": "B",
                    "concentrationStr": "1e-08",
                    "concentration": 1e-8,
                    "reads": 90,
                },
            ]
        )
        pooled = compute_frequency_signal(frame)
        apart = compute_frequency_signal(frame, REPS)
        assert pooled["signal"].to_list() == [pytest.approx(0.1), pytest.approx(0.9)]
        assert apart["signal"].to_list() == [1.0, 1.0]


class TestMeanBinFrame:
    def _signal(self) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "clonotypeKey": ["v1", "v1"],
                "replicate": ["B", "A"],
                "concentrationStr": ["1e-08", "1e-08"],
                "concentration": [1e-8, 1e-8],
                "signal": [2.0, 1.0],
            }
        )

    def test_without_replicate_keys_the_frame_keeps_its_two_key_columns(self):
        assert build_mean_bin_frame(self._signal()).columns == ["clonotypeKey", "concentrationStr", "meanBin"]

    def test_with_replicate_keys_the_replicate_joins_the_key_and_the_sort(self):
        frame = build_mean_bin_frame(self._signal(), REPS)
        assert frame.columns == ["clonotypeKey", "replicate", "concentrationStr", "meanBin"]
        assert frame["replicate"].to_list() == ["A", "B"]
