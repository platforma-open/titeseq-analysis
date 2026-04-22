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
    validate_sort_fraction,
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

    # Non-finite concentrations reach curve_fit as garbage inputs — downstream failures
    # are opaque and per-clonotype. Catching here is the only place with enough context
    # to produce an actionable error.
    def test_nan_concentration_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s",
                    "concentrationStr": "nan",
                    "concentration": float("nan"),
                    "bin": 1,
                    "reads": 5,
                },
            ]
        )
        with pytest.raises(InputValidationError, match="non-finite"):
            validate_concentration_column(df, has_bin=True)

    def test_inf_concentration_raises(self):
        df = _mk(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s",
                    "concentrationStr": "inf",
                    "concentration": float("inf"),
                    "bin": 1,
                    "reads": 5,
                },
            ]
        )
        with pytest.raises(InputValidationError, match="non-finite"):
            validate_concentration_column(df, has_bin=True)

    def test_concentration_over_ceiling_raises_actionable_error(self):
        # Above the int64 attomolar ceiling (~9.22 M) — any such value almost
        # certainly indicates a unit-entry mistake. Validation should flag it
        # with a clear error instead of letting it reach the attomolar cast.
        df = _mk(
            [
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
        with pytest.raises(InputValidationError, match="attomolar-encoding ceiling"):
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

    def test_all_empty_bin_raises_friendly_error(self):
        # An Integer metadata column whose sample values are all null arrives as
        # a String column of empty strings (TSV null encoding). The user should
        # see a pick-a-different-column message, not a raw "got String" dtype
        # error that doesn't hint at the root cause.
        df = pl.DataFrame(
            {
                "clonotypeKey": ["A", "A"],
                "sampleId": ["s1", "s2"],
                "concentrationStr": ["1", "1"],
                "concentration": [1.0, 1.0],
                "bin": ["", ""],
                "reads": [5, 5],
            },
            schema={
                "clonotypeKey": pl.Utf8,
                "sampleId": pl.Utf8,
                "concentrationStr": pl.Utf8,
                "concentration": pl.Float64,
                "bin": pl.Utf8,
                "reads": pl.Int64,
            },
        )
        with pytest.raises(InputValidationError, match="empty for every sample"):
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

    def test_all_empty_concentration_raises_friendly_error(self):
        # When the selected concentration metadata column is all-null, the TSV builder
        # emits empty strings and polars types the column as Utf8. Cast-to-Float64 would
        # surface as an opaque polars error — we raise a friendlier message first.
        df = pl.DataFrame(
            [
                {"clonotypeKey": "A", "sampleId": "s1", "concentration": "", "bin": 1, "reads": 5},
                {"clonotypeKey": "A", "sampleId": "s2", "concentration": "", "bin": 1, "reads": 5},
            ]
        )
        with pytest.raises(InputValidationError, match="no numeric values"):
            canonicalize_concentration(df)


class TestNarrowConcentrationRangeWarning:
    """R5 design-level guardrail: warn if non-zero concentrations span < 1 order of magnitude.

    A narrow dose range may not bracket Kd,app for any antibody, yielding
    kdOutOfRange = true on all clonotypes. This check does not block execution;
    it surfaces as a warning.
    """

    # 0 M (no-antigen control) must be excluded from the max/min ratio — otherwise
    # any dataset with a 0 M point would always be flagged.
    # Concentrations are in molar; values chosen to stay below the R2 ceiling
    # (~9.22 M) while preserving the max/min ratios the warning depends on.
    @pytest.mark.parametrize(
        "concs, expect_warning",
        [
            ([1e-9, 2e-9, 3e-9, 5e-9, 7e-9], True),    # 7/1 = 7 < 10 → warn
            ([1e-9, 1e-8, 1e-7], False),                # 100/1 = 100 ≥ 10 → no warn
            ([1e-9, 1e-8], False),                      # 10/1 = 10 exactly at boundary → no warn
            ([1e-9, 9.99e-9], True),                    # 9.99/1 = 9.99 < 10 → warn
            ([0.0, 1e-9, 1e-8, 1e-7], False),           # 0 excluded from ratio → 100/1 ok
            ([0.0, 1e-9, 2e-9, 3e-9], True),            # 0 excluded → 3/1 < 10 → warn
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


from io_layer import validate_bin_concentration_grid as _validate_grid


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


def _mk_reads_with_sort_fraction(
    concentrations: list[str],
    sort_fractions_by_conc: dict[str, list[float]],
    *,
    clonotypes: list[str] | None = None,
) -> pl.DataFrame:
    """Build a reads frame with per-sample sort_fraction metadata.

    One sample per (conc, bin); reads are duplicated across `clonotypes` so the
    validator's dedupe-by-sample logic has something non-trivial to chew on.
    """
    clonotypes = clonotypes or ["A"]
    rows = []
    for conc in concentrations:
        fractions = sort_fractions_by_conc[conc]
        for j, frac in enumerate(fractions):
            bin_label = j + 1
            for clone in clonotypes:
                rows.append(
                    {
                        "clonotypeKey": clone,
                        "sampleId": f"s_c{conc}_b{bin_label}",
                        "concentrationStr": conc,
                        "concentration": float(conc),
                        "bin": bin_label,
                        "reads": 10,
                        "sort_fraction": frac,
                    }
                )
    return pl.DataFrame(rows)


class TestValidateSortFraction:
    """FACS sort-fraction column validator — ingestion gate for the correction path."""

    # Happy path: sum-to-one within tolerance at every concentration, all values in-range.
    # The dedupe-by-sample is load-bearing here: passing the per-clonotype raw reads
    # frame would multiply every fraction by the number of clonotypes and falsely
    # inflate the per-concentration sum.
    def test_happy_path(self):
        df = _mk_reads_with_sort_fraction(
            concentrations=["1", "10"],
            sort_fractions_by_conc={
                "1": [0.25, 0.25, 0.25, 0.25],
                "10": [0.5, 0.1, 0.2, 0.2],
            },
            clonotypes=["A", "B", "C"],
        )
        validate_sort_fraction(df, "sort_fraction")

    def test_missing_column_raises(self):
        df = pl.DataFrame(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 5,
                }
            ]
        )
        with pytest.raises(InputValidationError, match="missing"):
            validate_sort_fraction(df, "sort_fraction")

    @pytest.mark.parametrize("bad_value", [-0.01, 1.01, 2.0, -5.0])
    def test_out_of_range_raises(self, bad_value):
        # In-bounds fractions for every other bin so only the one offender triggers.
        df = _mk_reads_with_sort_fraction(
            concentrations=["1"],
            sort_fractions_by_conc={"1": [bad_value, 1.0 - bad_value]},
        )
        with pytest.raises(InputValidationError, match=r"\[0, 1\]"):
            validate_sort_fraction(df, "sort_fraction")

    # User needs to know WHICH concentration is misconfigured — a generic
    # "sums don't match" message forces the user to rediscover the offender
    # themselves across potentially dozens of concentrations.
    def test_sum_violation_names_offending_concentration(self):
        df = _mk_reads_with_sort_fraction(
            concentrations=["1", "10"],
            sort_fractions_by_conc={
                "1": [0.25, 0.25, 0.25, 0.25],
                "10": [0.5, 0.3, 0.1, 0.0],  # sums to 0.9
            },
        )
        with pytest.raises(InputValidationError, match="10"):
            validate_sort_fraction(df, "sort_fraction")

    # Tolerance boundary: |sum − 1| < 1e-3 passes, ≥ 1e-3 fails. Confirms the
    # tolerance plane isn't off-by-one or mis-signed.
    def test_sum_tolerance_just_inside_passes(self):
        df = _mk_reads_with_sort_fraction(
            concentrations=["1"],
            sort_fractions_by_conc={"1": [0.2505, 0.2505, 0.2495, 0.2495]},  # sum=1.0000
        )
        validate_sort_fraction(df, "sort_fraction")

    def test_sum_tolerance_just_outside_raises(self):
        df = _mk_reads_with_sort_fraction(
            concentrations=["1"],
            sort_fractions_by_conc={"1": [0.252, 0.252, 0.250, 0.250]},  # sum=1.004 > 1.001
        )
        with pytest.raises(InputValidationError, match="sum to 1.0"):
            validate_sort_fraction(df, "sort_fraction")

    # A join miss produces null; the user must see it at validation, not lose
    # reads silently through polars' null-skipping group_by aggregation.
    def test_null_raises(self):
        df = pl.DataFrame(
            [
                {
                    "clonotypeKey": "A",
                    "sampleId": "s_c1_b1",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 1,
                    "reads": 10,
                    "sort_fraction": None,
                },
                {
                    "clonotypeKey": "A",
                    "sampleId": "s_c1_b2",
                    "concentrationStr": "1",
                    "concentration": 1.0,
                    "bin": 2,
                    "reads": 10,
                    "sort_fraction": 1.0,
                },
            ]
        )
        with pytest.raises(InputValidationError, match="null"):
            validate_sort_fraction(df, "sort_fraction")
