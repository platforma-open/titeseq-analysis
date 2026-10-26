"""Input/output layer (R1-R5): read reads table, validate schema, canonicalize concentrations."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from constants import (
    COL_ANTIGEN,
    COL_BIN,
    COL_CLONOTYPE,
    COL_CONC_STR,
    COL_CONC_VAL,
    COL_READS,
    COL_SAMPLE,
)


class InputValidationError(ValueError):
    """Raised when the input table violates R1-R5."""


def read_reads_table(path: str | Path) -> pl.DataFrame:
    """Load reads frame from parquet or tsv; dispatch by extension."""
    p = Path(path)
    if p.suffix in (".parquet", ".pq"):
        return pl.read_parquet(p)
    return pl.read_csv(p, separator="\t")


def validate_reads_schema(df: pl.DataFrame, has_bin: bool, has_antigen: bool) -> None:
    """R1, R3, R4: verify required columns exist for the active mode."""
    required = {COL_CLONOTYPE, COL_SAMPLE, COL_CONC_STR, COL_CONC_VAL, COL_READS}
    if has_bin:
        required.add(COL_BIN)
    if has_antigen:
        required.add(COL_ANTIGEN)
    missing = required - set(df.columns)
    if missing:
        raise InputValidationError(f"missing required columns: {sorted(missing)}")


def validate_concentration_column(df: pl.DataFrame, has_bin: bool) -> list[str]:
    """R2: concentrations must be finite non-negative floats. Unit-agnostic.

    Per spec R2 the block treats concentration values as dimensionless floats —
    no magnitude ceiling. The user's column may be in any unit (M, nM, µM, …).

    Returns list of warning strings (e.g. 0 M control without bin assignment in
    bin mode; R5 narrow-range warning when non-zero max/min spans less than one
    order of magnitude). Raises InputValidationError on any value < 0 or
    non-finite (NaN/Inf).
    """
    exprs = [
        (pl.col(COL_CONC_VAL) < 0).sum().alias("neg"),
        pl.col(COL_CONC_VAL).is_nan().sum().alias("nan"),
        pl.col(COL_CONC_VAL).is_infinite().sum().alias("inf"),
    ]
    if has_bin:
        exprs.append(((pl.col(COL_CONC_VAL) == 0) & pl.col(COL_BIN).is_null()).sum().alias("null_bin_at_zero"))
    counts = df.select(exprs).row(0, named=True)

    if counts["neg"] > 0:
        raise InputValidationError("concentration contains negative values")
    if counts["nan"] > 0 or counts["inf"] > 0:
        raise InputValidationError(
            f"concentration contains non-finite values (NaN: {counts['nan']}, Inf: {counts['inf']}); "
            "downstream curve fitting requires finite concentrations"
        )

    warnings: list[str] = []
    if has_bin and counts["null_bin_at_zero"] > 0:
        warnings.append(
            f"{counts['null_bin_at_zero']} rows at concentration 0 lack a bin assignment; "
            "ambiguous 0 M control entries (R2)"
        )

    # R5: narrow dose range may not bracket Kd,app; excludes 0 M control.
    nonzero = df.filter(pl.col(COL_CONC_VAL) > 0).select(pl.col(COL_CONC_VAL).unique())
    n_nonzero = nonzero.height
    if n_nonzero == 1:
        warnings.append("only one non-zero concentration present; narrow dose range may not bracket Kd,app (R5)")
    elif n_nonzero >= 2:
        max_c = float(nonzero.select(pl.col(COL_CONC_VAL).max()).item())
        min_c = float(nonzero.select(pl.col(COL_CONC_VAL).min()).item())
        if min_c > 0 and max_c / min_c < 10.0:
            warnings.append(
                f"non-zero concentrations span {max_c / min_c:.2f}x (less than one order of magnitude); "
                "narrow dose range may not bracket Kd,app (R5)"
            )
    return warnings


def _sorted_rows(reads: pl.DataFrame) -> pl.DataFrame:
    """Rows with a concentration and a bin: the only rows a sort correction applies to. A pre-sort
    library sample or an unsorted 0 M control may leave the per-sample cell metadata empty."""
    return reads.filter(pl.col(COL_CONC_VAL).is_not_null() & pl.col(COL_BIN).is_not_null())


def validate_sort_fraction(reads: pl.DataFrame, col: str) -> None:
    """FACS sort-fraction column validator (Adams, Mora, Walczak, Kinney 2016 eq. A3).

    Invariants:
      * Column exists in `reads`.
      * No nulls on sorted rows (a concentration and a bin); other rows are not corrected.
      * Values in [0, 1].
      * Per concentration, |Σ_b sort_fraction − 1| < 1e-3 — deduplicated by sample
        so that repeating a per-sample value across clonotypes does not inflate
        the per-concentration sum.

    No silent renormalization. Violations raise with a concrete message naming
    the offending concentration so the user knows which row to fix.
    """
    if col not in reads.columns:
        raise InputValidationError(f"Sort-fraction column '{col}' missing from reads")

    reads = _sorted_rows(reads)
    vals = reads[col]
    if vals.is_null().any():
        raise InputValidationError(f"Sort-fraction column '{col}' has null values")
    if (vals < 0).any() or (vals > 1).any():
        raise InputValidationError(f"Sort-fraction values must be in [0, 1]; column '{col}'")

    # `sort_fraction` is a per-sample metadata column joined across every clonotype
    # in that sample, so summing over the raw reads frame would count each sample
    # once per clonotype. Deduplicate by (sampleId, concentrationStr, col) first.
    per_sample = reads.select([COL_SAMPLE, COL_CONC_STR, col]).unique()
    sums = per_sample.group_by(COL_CONC_STR, maintain_order=True).agg(pl.col(col).sum().alias("s"))
    bad = sums.filter((pl.col("s") - 1.0).abs() > 1e-3)
    if bad.height > 0:
        violations = bad.select([COL_CONC_STR, "s"]).to_dicts()
        raise InputValidationError(
            f"sort_fraction values must sum to 1.0 per concentration (tolerance 1e-3). Violations: {violations}"
        )


def validate_cell_count(reads: pl.DataFrame, col: str) -> None:
    """Cells sorted into each sample (Adams, Mora, Walczak, Kinney 2016 eq. E4, C_bc).

    Invariants on sorted rows (a concentration and a bin): no nulls, finite, >= 0, one value per
    sample, and at least one sample with cells at every concentration. The column must be present.
    """
    if col not in reads.columns:
        raise InputValidationError(f"Cell-count column '{col}' missing from reads")
    reads = _sorted_rows(reads)
    vals = reads[col].cast(pl.Float64)
    if vals.is_null().any():
        raise InputValidationError(f"Cell-count column '{col}' has null values")
    if (~vals.is_finite()).any() or (vals < 0).any():
        raise InputValidationError(f"Cell counts must be finite and >= 0; column '{col}'")
    per_sample = reads.select([COL_SAMPLE, col]).unique()
    dup = per_sample.group_by(COL_SAMPLE).len().filter(pl.col("len") > 1)
    if dup.height > 0:
        raise InputValidationError(
            f"Cell count must be one value per sample; samples with several: {dup[COL_SAMPLE].to_list()[:5]}"
        )
    per_conc = reads.select([COL_SAMPLE, COL_CONC_STR, col]).unique().group_by(COL_CONC_STR).agg(pl.col(col).sum())
    empty = per_conc.filter(pl.col(col) <= 0)
    if empty.height > 0:
        raise InputValidationError(f"No cells sorted at concentration(s) {empty[COL_CONC_STR].to_list()}")


def validate_bin_column(df: pl.DataFrame) -> None:
    """R3: bin must be a positive integer. Non-consecutive labels are allowed.

    Nulls are permitted at concentration 0 (warned separately); `<=0` on null
    yields null which is ignored by `.sum()`, so no intermediate filter needed.
    """
    # An Integer metadata column whose sample values are all null arrives here as a
    # string column of empty strings (the TSV builder's null encoding). Surface a
    # user-actionable message instead of the generic "got String" dtype complaint.
    if df.schema[COL_BIN] == pl.Utf8:
        non_empty = df.filter(pl.col(COL_BIN).is_not_null() & (pl.col(COL_BIN).str.len_chars() > 0)).height
        if non_empty == 0:
            raise InputValidationError(
                "bin column has no integer values — the selected metadata "
                "column appears to be empty for every sample. Pick a bin "
                "column that is populated across samples (e.g. a numeric "
                "FACS bin-label column)."
            )
    if not df[COL_BIN].dtype.is_integer():
        raise InputValidationError(f"bin column must be integer type, got {df[COL_BIN].dtype}")
    if df.select((pl.col(COL_BIN) <= 0).sum()).item() > 0:
        raise InputValidationError("bin column must contain positive integers (>= 1)")


def max_bin_label(df: pl.DataFrame) -> int:
    """Return the actual max bin label present. Spec: label, not count — labels may be non-consecutive."""
    return int(df.select(pl.col(COL_BIN).max()).item())


def validate_antigen_filter(df: pl.DataFrame, antigen_column_ref: str | None, target_antigen: str | None) -> list[str]:
    """R4: antigen filtering semantics.

    - antigenColumnRef without targetAntigen -> user-facing error.
    - targetAntigen without antigenColumnRef -> warning (single-antigen assumption).
    - targetAntigen not present in antigen column -> error.
    """
    warnings: list[str] = []
    if antigen_column_ref is not None and target_antigen is None:
        raise InputValidationError("antigenColumnRef provided without targetAntigen — specify which antigen to analyse")
    if antigen_column_ref is None and target_antigen is not None:
        warnings.append("targetAntigen set without antigenColumnRef — treating dataset as single-antigen")
        return warnings
    if antigen_column_ref is not None and target_antigen is not None:
        present = set(df[COL_ANTIGEN].drop_nulls().unique().to_list())
        if target_antigen not in present:
            raise InputValidationError(
                f"targetAntigen {target_antigen!r} not found in {COL_ANTIGEN} column; present values: {sorted(present)}"
            )
    return warnings


def apply_antigen_filter(df: pl.DataFrame, target_antigen: str | None) -> pl.DataFrame:
    """R4: keep only rows matching target_antigen (if antigen column present)."""
    if target_antigen is None or COL_ANTIGEN not in df.columns:
        return df
    return df.filter(pl.col(COL_ANTIGEN) == target_antigen)


def validate_sample_metadata_uniqueness(df: pl.DataFrame, has_bin: bool, has_antigen: bool) -> None:
    """R5: each sampleId must map to a unique (concentration, bin, antigen) tuple.

    Metadata is constant for every clonotype in a sample, so a full-table group_by
    hashes millions of rows to recover ≤ n_samples unique tuples. Dedupe the key
    columns first — this keeps the check O(samples) instead of O(rows).
    """
    key_cols = [COL_CONC_STR]
    if has_bin:
        key_cols.append(COL_BIN)
    if has_antigen:
        key_cols.append(COL_ANTIGEN)

    sample_metadata = df.select([COL_SAMPLE, *key_cols]).unique()
    per_sample_unique = sample_metadata.group_by(COL_SAMPLE, maintain_order=True).agg(
        [pl.col(c).n_unique().alias(f"{c}_nunique") for c in key_cols]
    )
    for c in key_cols:
        offenders = per_sample_unique.filter(pl.col(f"{c}_nunique") > 1)
        if offenders.height > 0:
            sample_list = offenders[COL_SAMPLE].to_list()[:5]
            raise InputValidationError(f"sampleId must map to a single {c}; violating samples: {sample_list}")


def validate_bin_concentration_grid(samples: pl.DataFrame | None) -> list[str]:
    """R5 sub-clause: warn when a (bin, concentration) sample is absent from the titration.

    An unmeasured (bin, concentration) sample drops a bin from mean_bin at that
    concentration, which biases the result with no other flag.

    `samples` is the sample table — one row per sampleId with its bin and
    concentration — not the reads table. The reads table cannot answer this
    question. It is sparse, so a (bin, concentration) pair is missing from it both
    when the sample was never measured and when every clonotype simply had no reads
    there. The second case is routine and unbiased: at 0 M the high-fluorescence
    bins are empty by construction, and a clonotype absent from a bin is a true zero
    that contributes nothing to either side of the mean_bin ratio. Checking reads
    therefore fires on every real dataset (bug-hunt M3). `None` skips the check.
    """
    warnings: list[str] = []
    if samples is None or COL_BIN not in samples.columns:
        return warnings

    measured = samples.select([COL_BIN, COL_CONC_STR]).unique()
    bins = samples[COL_BIN].unique()
    concentrations = samples[COL_CONC_STR].unique()
    expected = len(bins) * len(concentrations)
    if expected <= 1 or measured.height >= expected:
        return warnings

    absent = (
        bins.to_frame(name=COL_BIN)
        .join(concentrations.to_frame(name=COL_CONC_STR), how="cross")
        .join(measured, on=[COL_BIN, COL_CONC_STR], how="anti")
        .sort([COL_CONC_STR, COL_BIN])
    )
    examples = [f"(bin {b}, {c} M)" for b, c in absent.head(5).iter_rows()]
    warnings.append(
        f"{absent.height} (bin, concentration) sample(s) absent from the titration; "
        f"expected {expected}, examples: {examples} (R5)"
    )
    return warnings


def canonicalize_concentration(df: pl.DataFrame) -> pl.DataFrame:
    """Ensure both `concentrationStr` (canonical) and `concentration` (Float64) present (R14).

    `concentrationStr` carries the byte-exact upstream metadata value; it is what
    becomes the `pl7.app/vdj/concentration` String axis in output PColumns. The
    Tengo workflow's xsv.importFile wraps this column directly as a String axis
    (no parse step). `concentration` is the Float64 we use internally for fit
    arithmetic — never written to output, to avoid a float→string→float
    roundtrip between Python and the workflow.
    """
    if COL_CONC_STR not in df.columns:
        df = df.with_columns(pl.col(COL_CONC_VAL).cast(pl.Utf8).alias(COL_CONC_STR))
    if COL_CONC_VAL not in df.columns:
        df = df.with_columns(pl.col(COL_CONC_STR).cast(pl.Float64).alias(COL_CONC_VAL))
    if df.schema[COL_CONC_VAL] != pl.Float64:
        # Same all-null-column pattern handled in `validate_bin_column`.
        if df.schema[COL_CONC_VAL] == pl.Utf8:
            non_empty = df.filter(
                pl.col(COL_CONC_VAL).is_not_null() & (pl.col(COL_CONC_VAL).str.len_chars() > 0)
            ).height
            if non_empty == 0:
                raise InputValidationError(
                    "concentration column has no numeric values — the selected metadata "
                    "column appears to be empty for every sample. Pick a concentration "
                    "column that is populated across samples (e.g. a numeric molar column)."
                )
        df = df.with_columns(pl.col(COL_CONC_VAL).cast(pl.Float64))
    return df
