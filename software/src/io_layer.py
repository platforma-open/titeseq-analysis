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
    MAX_CONCENTRATION_M,
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
    """R2: concentrations must be non-negative floats within the attomolar-encodable range.

    Returns list of warning strings (e.g. 0 M control without bin assignment in bin mode;
    R5 narrow-range warning when non-zero max/min spans less than one order of magnitude).
    Raises InputValidationError on any value < 0 or above MAX_CONCENTRATION_M (int64
    ceiling for attomolar encoding; see constants.py).
    """
    exprs = [
        (pl.col(COL_CONC_VAL) < 0).sum().alias("neg"),
        (pl.col(COL_CONC_VAL) > MAX_CONCENTRATION_M).sum().alias("over_max"),
    ]
    if has_bin:
        exprs.append(((pl.col(COL_CONC_VAL) == 0) & pl.col(COL_BIN).is_null()).sum().alias("null_bin_at_zero"))
    counts = df.select(exprs).row(0, named=True)

    if counts["neg"] > 0:
        raise InputValidationError("concentration contains negative values")
    if counts["over_max"] > 0:
        offenders = (
            df.filter(pl.col(COL_CONC_VAL) > MAX_CONCENTRATION_M)[COL_CONC_VAL].unique().sort(descending=True).to_list()
        )
        raise InputValidationError(
            f"concentration exceeds {MAX_CONCENTRATION_M:.3g} M (attomolar-encoding ceiling); "
            f"offending value(s): {offenders[:5]}. "
            "TiteSeq concentrations are typically sub-µM — check the unit column and metadata entry."
        )

    warnings: list[str] = []
    if has_bin and counts["null_bin_at_zero"] > 0:
        warnings.append(
            f"{counts['null_bin_at_zero']} rows at concentration 0 lack a bin assignment; "
            "ambiguous 0 M control entries (R2)"
        )

    # R5: narrow dose range may not bracket K_D,app; excludes 0 M control.
    nonzero = df.filter(pl.col(COL_CONC_VAL) > 0).select(pl.col(COL_CONC_VAL).unique())
    n_nonzero = nonzero.height
    if n_nonzero == 1:
        warnings.append(
            "only one non-zero concentration present; narrow dose range may not bracket K_D,app (R5)"
        )
    elif n_nonzero >= 2:
        max_c = float(nonzero.select(pl.col(COL_CONC_VAL).max()).item())
        min_c = float(nonzero.select(pl.col(COL_CONC_VAL).min()).item())
        if min_c > 0 and max_c / min_c < 10.0:
            warnings.append(
                f"non-zero concentrations span {max_c / min_c:.2f}x (less than one order of magnitude); "
                "narrow dose range may not bracket K_D,app (R5)"
            )
    return warnings


def validate_bin_column(df: pl.DataFrame) -> None:
    """R3: bin must be a positive integer. Non-consecutive labels are allowed.

    Nulls are permitted at concentration 0 (warned separately); `<=0` on null
    yields null which is ignored by `.sum()`, so no intermediate filter needed.
    """
    if df[COL_BIN].dtype not in (pl.Int8, pl.Int16, pl.Int32, pl.Int64, pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64):
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
    """R5: each sampleId must map to a unique (concentration, bin, antigen) tuple."""
    key_cols = [COL_CONC_STR]
    if has_bin:
        key_cols.append(COL_BIN)
    if has_antigen:
        key_cols.append(COL_ANTIGEN)

    per_sample_unique = df.group_by(COL_SAMPLE).agg([pl.col(c).n_unique().alias(f"{c}_nunique") for c in key_cols])
    for c in key_cols:
        offenders = per_sample_unique.filter(pl.col(f"{c}_nunique") > 1)
        if offenders.height > 0:
            sample_list = offenders[COL_SAMPLE].to_list()[:5]
            raise InputValidationError(f"sampleId must map to a single {c}; violating samples: {sample_list}")


def validate_bin_concentration_grid(df: pl.DataFrame) -> list[str]:
    """R5 sub-clause: warn when (bin, concentration) combos are non-uniformly populated.

    Absent combinations silently reduce the number of bins contributing to mean_bin
    at that concentration, biasing the result without any flag. The reference grid
    is every (bin, concentrationStr) pair observed for ANY clonotype; any clonotype
    with fewer combinations present is flagged.
    """
    warnings: list[str] = []
    if COL_BIN not in df.columns:
        return warnings

    ref_grid = df.select([COL_BIN, COL_CONC_STR]).unique()
    n_grid = ref_grid.height
    if n_grid <= 1:
        return warnings

    per_clone = (
        df.select([COL_CLONOTYPE, COL_BIN, COL_CONC_STR])
        .unique()
        .group_by(COL_CLONOTYPE)
        .agg(pl.len().alias("n_combos"))
    )
    missing = per_clone.filter(pl.col("n_combos") < n_grid)
    if missing.height > 0:
        offenders = missing[COL_CLONOTYPE].to_list()[:5]
        warnings.append(
            f"{missing.height} clonotype(s) have missing (bin, concentration) combinations; "
            f"expected {n_grid}, examples: {offenders} (R5)"
        )
    return warnings


def canonicalize_concentration(df: pl.DataFrame) -> pl.DataFrame:
    """Ensure both concentration axes present; preserve original string exactly (R14).

    `concentrationStr` is the canonical internal join key — it must compare equal between
    rows that came from the same upstream metadata value, even when float parsing would
    introduce drift. Preserving the input string avoids float→string→float roundtrips.
    The output PColumn axis is `concentrationAM` (Long), so `concentrationStr` is no longer
    surfaced to Graph Maker and its lexicographic sortability is irrelevant.
    """
    if COL_CONC_STR not in df.columns:
        df = df.with_columns(pl.col(COL_CONC_VAL).cast(pl.Utf8).alias(COL_CONC_STR))
    if COL_CONC_VAL not in df.columns:
        df = df.with_columns(pl.col(COL_CONC_STR).cast(pl.Float64).alias(COL_CONC_VAL))
    if df.schema[COL_CONC_VAL] != pl.Float64:
        # A Float metadata column whose sample values are all null arrives here as a
        # string column of empty strings (the TSV builder's null encoding). Detect
        # this and raise a clear error instead of a raw polars cast traceback.
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
