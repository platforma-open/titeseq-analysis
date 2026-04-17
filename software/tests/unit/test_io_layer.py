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
        df = _mk(
            [
                {"sampleId": "s1", "concentrationStr": "1", "concentration": 1.0, "bin": 1, "reads": 5},
            ]
        )
        with pytest.raises(InputValidationError, match="clonotypeKey"):
            validate_reads_schema(df, has_bin=True, has_antigen=False)

    def test_missing_sample_col_raises(self):
        df = _mk(
            [
                {"clonotypeKey": "A", "concentrationStr": "1", "concentration": 1.0, "bin": 1, "reads": 5},
            ]
        )
        with pytest.raises(InputValidationError, match="sampleId"):
            validate_reads_schema(df, has_bin=True, has_antigen=False)


class TestConcentrationValidation:
    def test_negative_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s",
                    "concentrationStr": "-1",
                    "concentration": -1.0,
                    "bin": 1,
                    "reads": 5,
                },
            ]
        )
        with pytest.raises(InputValidationError, match="negative"):
            validate_concentration_column(df, has_bin=True)


class TestBinValidation:
    def test_non_consecutive_labels_ok(self):
        # Spec: labels may be non-consecutive; max_bin_label is actual max, not count.
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": b,
                    "reads": 5,
                }
                for b in [1, 2, 5, 8]
            ]
        )
        validate_bin_column(df)
        assert max_bin_label(df) == 8  # label, not count

    def test_zero_bin_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 0,
                    "reads": 5,
                },
            ]
        )
        with pytest.raises(InputValidationError, match="positive"):
            validate_bin_column(df)

    def test_negative_bin_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": -1,
                    "reads": 5,
                },
            ]
        )
        with pytest.raises(InputValidationError, match="positive"):
            validate_bin_column(df)

    def test_non_integer_bin_raises(self):
        df = pl.DataFrame(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 1.5,
                    "reads": 5,
                }
            ],
        )
        with pytest.raises(InputValidationError, match="integer"):
            validate_bin_column(df)


class TestAntigenFilter:
    def test_ref_without_target_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 5,
                    "antigen": "X",
                }
            ]
        )
        with pytest.raises(InputValidationError, match="targetAntigen"):
            validate_antigen_filter(df, antigen_column_ref="antigen", target_antigen=None)

    def test_target_without_ref_warns(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 5,
                }
            ]
        )
        warnings = validate_antigen_filter(df, antigen_column_ref=None, target_antigen="X")
        assert any("antigenColumnRef" in w for w in warnings)

    def test_target_not_in_column_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 5,
                    "antigen": "Y",
                }
            ]
        )
        with pytest.raises(InputValidationError, match="not found"):
            validate_antigen_filter(df, antigen_column_ref="antigen", target_antigen="X")


class TestSampleUniqueness:
    def test_sample_with_two_concentrations_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 5,
                },
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "10",
                    "concentration": 10.0,
                    "bin": 1,
                    "reads": 5,
                },
            ]
        )
        with pytest.raises(InputValidationError, match="concentrationStr"):
            validate_sample_metadata_uniqueness(df, has_bin=True, has_antigen=False)

    def test_sample_with_two_bins_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 5,
                },
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 2,
                    "reads": 5,
                },
            ]
        )
        with pytest.raises(InputValidationError, match="bin"):
            validate_sample_metadata_uniqueness(df, has_bin=True, has_antigen=False)


class TestCanonicalConcentrationAxis:
    # R14: string-preserving canonicalization prevents float-drift in axis joins.
    def test_numeric_equal_but_string_different_kept_separate(self):
        df = pl.DataFrame(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "1.0",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 5,
                },
                {
                    "clonotypeKey": "A",
                    "sampleId": "s2",
                    "concentrationStr": "1.000",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 5,
                },
            ]
        )
        result = canonicalize_concentration(df)
        assert set(result["concentrationStr"].to_list()) == {"1.0", "1.000"}

    def test_adds_missing_concentrationStr(self):
        df = pl.DataFrame([{"clonotypeKey": "A", "sampleId": "s", "concentration": 0.001, "bin": 1, "reads": 5}])
        result = canonicalize_concentration(df)
        assert "concentrationStr" in result.columns
        assert result["concentrationStr"][0] == "0.001"


class TestNarrowConcentrationRangeWarning:
    """R5 design-level guardrail: warn if non-zero concentrations span < 1 order of magnitude.

    A narrow dose range may not bracket K_D,app for any antibody, yielding
    kdOutOfRange = true on all clonotypes. This check does not block execution;
    it surfaces as a warning.
    """

    # 0 M (no-antigen control) must be excluded from the max/min ratio — otherwise
    # any dataset with a 0 M point would always be flagged.
    @pytest.mark.parametrize(
        "concs, expect_warning",
        [
            ([1.0, 2.0, 3.0, 5.0, 7.0], True),    # 7/1 = 7 < 10 → warn
            ([1.0, 10.0, 100.0], False),          # 100/1 = 100 ≥ 10 → no warn
            ([1.0, 10.0], False),                 # 10/1 = 10 exactly at boundary → no warn
            ([1.0, 9.99], True),                  # 9.99/1 = 9.99 < 10 → warn
            ([0.0, 1.0, 10.0, 100.0], False),     # 0 excluded from ratio → 100/1 ok
            ([0.0, 1.0, 2.0, 3.0], True),         # 0 excluded → 3/1 < 10 → warn
        ],
    )
    def test_narrow_range_emits_warning(self, concs, expect_warning):
        rows = [
            {
                "clonotypeKey": "A",
                "sampleId": f"s{i}",
                "concentrationStr": str(c),
                "concentration": c,
                "bin": 1,
                "reads": 5,
            }
            for i, c in enumerate(concs)
        ]
        df = _mk(rows)
        warnings = validate_concentration_column(df, has_bin=True)
        has_narrow = any(
            ("order of magnitude" in w.lower()) or ("narrow" in w.lower()) or ("span" in w.lower())
            for w in warnings
        )
        assert has_narrow is expect_warning, f"concs={concs}, warnings={warnings}"

    # Single non-zero concentration has vacuously zero span → should warn.
    def test_single_nonzero_concentration_warns(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s1",
                    "concentrationStr": "0",
                    "concentration": 0.0,
                    "bin": 1,
                    "reads": 5,
                },
                {
                    "clonotypeKey": "A",
                    "sampleId": "s2",
                    "concentrationStr": "5",
                    "concentration": 5.0,
                    "bin": 1,
                    "reads": 5,
                },
            ]
        )
        warnings = validate_concentration_column(df, has_bin=True)
        assert any(
            ("narrow" in w.lower()) or ("range" in w.lower()) or ("span" in w.lower()) for w in warnings
        ), f"warnings={warnings}"


# R5 sub-clause: validate_bin_concentration_grid is not yet implemented.
# Once added to io_layer, these tests exercise it. Until then they are skipped.
try:
    from io_layer import validate_bin_concentration_grid as _validate_grid  # type: ignore

    _HAS_GRID_VALIDATOR = True
except ImportError:
    _HAS_GRID_VALIDATOR = False


@pytest.mark.skipif(
    not _HAS_GRID_VALIDATOR,
    reason="validate_bin_concentration_grid not implemented yet (R5 sub-clause)",
)
class TestBinConcentrationGrid:
    """R5 sub-clause: warn when (bin, concentration) combos are non-uniformly populated.

    Absent combinations silently reduce the number of bins contributing to
    mean_bin at that concentration, biasing the result without any flag.
    """

    def test_uniform_grid_no_warning(self):
        # Every clonotype has the full (bin × conc) cartesian product present.
        rows = []
        for clone in ["A", "B"]:
            for conc in ["1", "10"]:
                for b in [1, 2]:
                    rows.append(
                        {
                            "clonotypeKey": clone,
                            "sampleId": f"s-{clone}-{conc}-{b}",
                            "concentrationStr": conc,
                            "concentration": float(conc),
                            "bin": b,
                            "reads": 5,
                        }
                    )
        warnings = _validate_grid(_mk(rows))
        assert warnings == []

    def test_missing_combo_emits_warning(self):
        # A has (bin=1, conc=10); B is missing that combo → warn on B.
        rows = [
            {
                "clonotypeKey": "A",
                "sampleId": "s1",
                "concentrationStr": "1",
                "concentration": 1.0,
                "bin": 1,
                "reads": 5,
            },
            {
                "clonotypeKey": "A",
                "sampleId": "s2",
                "concentrationStr": "10",
                "concentration": 10.0,
                "bin": 1,
                "reads": 5,
            },
            {
                "clonotypeKey": "B",
                "sampleId": "s1",
                "concentrationStr": "1",
                "concentration": 1.0,
                "bin": 1,
                "reads": 5,
            },
        ]
        warnings = _validate_grid(_mk(rows))
        assert len(warnings) > 0
