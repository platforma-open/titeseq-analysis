"""Input/output layer (R1-R5): read reads table, validate schema, canonicalize concentrations."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from constants import (
    COL_ANTIGEN,
    COL_BIN,
    COL_BIN_FLUORESCENCE,
    COL_CLONOTYPE,
    COL_CONC_STR,
    COL_CONC_VAL,
    COL_READS,
    COL_SAMPLE,
    MAX_LOG10_FLUORESCENCE,
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


def validate_replicate_column(reads: pl.DataFrame, col: str) -> list[str]:
    """Replicate label, one value per sample. Optional input.

    Returns warnings rather than raising: unbound or unusable, the run proceeds as a single
    titration, which is what every project did before the column existed.
    """
    if col not in reads.columns:
        return [f"Replicate column '{col}' missing from reads; treating the run as one replicate"]
    labels = reads.get_column(col).cast(pl.Utf8, strict=False)
    if labels.is_null().any() or (labels.str.len_chars() == 0).any():
        return [f"Replicate column '{col}' has empty values; treating the run as one replicate"]
    per_sample = reads.select([COL_SAMPLE, col]).unique()
    dup = per_sample.group_by(COL_SAMPLE).len().filter(pl.col("len") > 1)
    if dup.height > 0:
        return [
            f"Replicate must be one value per sample; samples with several: {dup[COL_SAMPLE].to_list()[:5]}. "
            "Treating the run as one replicate"
        ]
    return []


def validate_sort_fraction(reads: pl.DataFrame, col: str, group_keys: list[str] | None = None) -> None:
    """FACS sort-fraction column validator (Adams, Mora, Walczak, Kinney 2016 eq. A3).

    Invariants:
      * Column exists in `reads`.
      * No nulls (every sample with reads must carry a fraction).
      * Values in [0, 1].
      * Per group, |Σ_b sort_fraction − 1| < 1e-3 — deduplicated by sample so that repeating a
        per-sample value across clonotypes does not inflate the per-group sum.

    `group_keys` names the columns the sum is taken within, concentration last. Each replicate
    sorts its own cells, so with replicates the group is (replicate, concentration): pooled, the
    fractions would sum to the replicate count.

    No silent renormalization. Violations raise with a concrete message naming
    the offending group so the user knows which row to fix.
    """
    if col not in reads.columns:
        raise InputValidationError(f"Sort-fraction column '{col}' missing from reads")

    vals = reads[col]
    if vals.is_null().any():
        raise InputValidationError(f"Sort-fraction column '{col}' has null values")
    if (vals < 0).any() or (vals > 1).any():
        raise InputValidationError(f"Sort-fraction values must be in [0, 1]; column '{col}'")

    keys = [COL_CONC_STR] if group_keys is None else group_keys
    # `sort_fraction` is a per-sample metadata column joined across every clonotype
    # in that sample, so summing over the raw reads frame would count each sample
    # once per clonotype. Deduplicate by (sampleId, group keys, col) first.
    per_sample = reads.select([COL_SAMPLE, *keys, col]).unique()
    sums = per_sample.group_by(keys, maintain_order=True).agg(pl.col(col).sum().alias("s"))
    bad = sums.filter((pl.col("s") - 1.0).abs() > 1e-3)
    if bad.height > 0:
        violations = bad.select([*keys, "s"]).to_dicts()
        scope = "concentration" if keys == [COL_CONC_STR] else " and ".join(keys)
        raise InputValidationError(
            f"sort_fraction values must sum to 1.0 per {scope} (tolerance 1e-3). Violations: {violations}"
        )


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


def validate_bin_fluorescence(df: pl.DataFrame) -> None:
    """Every row in a gate (bin >= 1) must carry a finite log10 gate position.

    A missing position would drop that gate's reads from the mean log10 fluorescence and move the
    signal without any warning, so it is refused. Rows with a null bin (allowed at concentration 0)
    are not checked.

    A position above MAX_LOG10_FLUORESCENCE is linear fluorescence, not log10. Fitted as log10 it
    moves K_D by decades, and above about 308 it overflows a float in the fit.
    """
    position = pl.col(COL_BIN_FLUORESCENCE).cast(pl.Float64, strict=False)
    gated = df.filter(pl.col(COL_BIN) >= 1)
    bad = gated.filter(position.is_null() | ~position.is_finite())
    if bad.height > 0:
        examples = bad.select([COL_SAMPLE, COL_CONC_STR, COL_BIN]).unique().sort(COL_SAMPLE).head(5).to_dicts()
        raise InputValidationError(
            f"Gate fluorescence column ('{COL_BIN_FLUORESCENCE}') must hold a finite log10 fluorescence "
            f"for every gate; {bad.height} rows are null or not finite, e.g. {examples}"
        )
    high = gated.filter(position > MAX_LOG10_FLUORESCENCE)
    if high.height > 0:
        top = high.select(position.max()).item()
        examples = high.select([COL_SAMPLE, COL_CONC_STR, COL_BIN]).unique().sort(COL_SAMPLE).head(5).to_dicts()
        raise InputValidationError(
            f"Gate fluorescence column ('{COL_BIN_FLUORESCENCE}') must hold log10 fluorescence, but it "
            f"reaches {top:g}, above the {MAX_LOG10_FLUORESCENCE:g} accepted. It looks like linear "
            f"fluorescence: supply log10 of each gate's mean. {high.height} rows, e.g. {examples}"
        )


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


def validate_bin_concentration_grid(df: pl.DataFrame) -> list[str]:
    """R5 sub-clause: report how many clonotypes miss part of the (bin, concentration) grid.

    Absent combinations reduce the number of bins contributing to mean_bin at that
    concentration. The reference grid is every (bin, concentrationStr) pair observed for ANY
    clonotype. Run it on the frame the block fits, after the antigen filter.

    The line carries a count and a share, no clonotype keys: keys are internal hashes, and a
    mutational-scanning library flags nearly every clonotype (1,150,360 on Adams 2016).
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
        .group_by(COL_CLONOTYPE, maintain_order=True)
        .agg(pl.len().alias("n_combos"))
    )
    missing = per_clone.filter(pl.col("n_combos") < n_grid).height
    if missing > 0:
        share = 100.0 * missing / per_clone.height
        warnings.append(
            f"{missing} of {per_clone.height} clonotypes ({share:.1f}%) are absent from at least one of "
            f"the {n_grid} (bin, concentration) samples (R5)"
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
