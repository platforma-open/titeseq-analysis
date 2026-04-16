"""Behavioral tests for io_layer.py (R1-R5 validation + antigen filter + canonicalization)."""

from __future__ import annotations

import polars as pl
import pytest

from io_layer import (
    InputValidationError,
    canonicalize_concentration,
    max_bin_label,
    validate_antigen_filter,
    validate_bin_column,
    validate_concentration_column,
    validate_reads_schema,
    validate_sample_metadata_uniqueness,
)


def _mk(rows, has_bin=True, has_antigen=False):
    return pl.DataFrame(rows)


class TestSchema:
    def test_missing_clonotype_col_raises(self):
        df = _mk([
            {"sampleId": "s1", "concentrationStr": "1",
             "concentration": 1.0, "bin": 1, "reads": 5},
        ])
        with pytest.raises(InputValidationError, match="clonotypeKey"):
            validate_reads_schema(df, has_bin=True, has_antigen=False)

    def test_missing_sample_col_raises(self):
        df = _mk([
            {"clonotypeKey": "A", "concentrationStr": "1",
             "concentration": 1.0, "bin": 1, "reads": 5},
        ])
        with pytest.raises(InputValidationError, match="sampleId"):
            validate_reads_schema(df, has_bin=True, has_antigen=False)


class TestConcentrationValidation:
    def test_negative_raises(self):
        df = _mk([
            {"clonotypeKey": "A", "sampleId": "s", "concentrationStr": "-1",
             "concentration": -1.0, "bin": 1, "reads": 5},
        ])
        with pytest.raises(InputValidationError, match="negative"):
            validate_concentration_column(df, has_bin=True)


class TestBinValidation:
    def test_non_consecutive_labels_ok(self):
        # Spec: labels may be non-consecutive; max_bin_label is actual max, not count.
        df = _mk([
            {"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1",
             "concentration": 1.0, "bin": b, "reads": 5}
            for b in [1, 2, 5, 8]
        ])
        validate_bin_column(df)
        assert max_bin_label(df) == 8  # label, not count

    def test_zero_bin_raises(self):
        df = _mk([
            {"clonotypeKey": "A", "sampleId": "s", "concentrationStr": "1",
             "concentration": 1.0, "bin": 0, "reads": 5},
        ])
        with pytest.raises(InputValidationError, match="positive"):
            validate_bin_column(df)

    def test_negative_bin_raises(self):
        df = _mk([
            {"clonotypeKey": "A", "sampleId": "s", "concentrationStr": "1",
             "concentration": 1.0, "bin": -1, "reads": 5},
        ])
        with pytest.raises(InputValidationError, match="positive"):
            validate_bin_column(df)

    def test_non_integer_bin_raises(self):
        df = pl.DataFrame(
            [{"clonotypeKey": "A", "sampleId": "s", "concentrationStr": "1",
              "concentration": 1.0, "bin": 1.5, "reads": 5}],
        )
        with pytest.raises(InputValidationError, match="integer"):
            validate_bin_column(df)


class TestAntigenFilter:
    def test_ref_without_target_raises(self):
        df = _mk([{"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1",
                   "concentration": 1.0, "bin": 1, "reads": 5, "antigen": "X"}])
        with pytest.raises(InputValidationError, match="targetAntigen"):
            validate_antigen_filter(df, antigen_column_ref="antigen", target_antigen=None)

    def test_target_without_ref_warns(self):
        df = _mk([{"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1",
                   "concentration": 1.0, "bin": 1, "reads": 5}])
        warnings = validate_antigen_filter(df, antigen_column_ref=None, target_antigen="X")
        assert any("antigenColumnRef" in w for w in warnings)

    def test_target_not_in_column_raises(self):
        df = _mk([{"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1",
                   "concentration": 1.0, "bin": 1, "reads": 5, "antigen": "Y"}])
        with pytest.raises(InputValidationError, match="not found"):
            validate_antigen_filter(df, antigen_column_ref="antigen", target_antigen="X")


class TestSampleUniqueness:
    def test_sample_with_two_concentrations_raises(self):
        df = _mk([
            {"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1",
             "concentration": 1.0, "bin": 1, "reads": 5},
            {"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "10",
             "concentration": 10.0, "bin": 1, "reads": 5},
        ])
        with pytest.raises(InputValidationError, match="concentrationStr"):
            validate_sample_metadata_uniqueness(df, has_bin=True, has_antigen=False)

    def test_sample_with_two_bins_raises(self):
        df = _mk([
            {"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1",
             "concentration": 1.0, "bin": 1, "reads": 5},
            {"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1",
             "concentration": 1.0, "bin": 2, "reads": 5},
        ])
        with pytest.raises(InputValidationError, match="bin"):
            validate_sample_metadata_uniqueness(df, has_bin=True, has_antigen=False)


class TestCanonicalConcentrationAxis:
    # R14: string-preserving canonicalization prevents float-drift in axis joins.
    def test_numeric_equal_but_string_different_kept_separate(self):
        df = pl.DataFrame(
            [
                {"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1.0",
                 "concentration": 1.0, "bin": 1, "reads": 5},
                {"clonotypeKey": "A", "sampleId": "s2", "concentrationStr": "1.000",
                 "concentration": 1.0, "bin": 1, "reads": 5},
            ]
        )
        result = canonicalize_concentration(df)
        assert set(result["concentrationStr"].to_list()) == {"1.0", "1.000"}

    def test_adds_missing_concentrationStr(self):
        df = pl.DataFrame(
            [{"clonotypeKey": "A", "sampleId": "s", "concentration": 0.001,
              "bin": 1, "reads": 5}]
        )
        result = canonicalize_concentration(df)
        assert "concentrationStr" in result.columns
        assert result["concentrationStr"][0] == "0.001"
