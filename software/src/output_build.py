"""R13/R14/R14b: assemble per-clonotype and per-(clonotype, concentration) output frames.

All axis joins are keyed on canonical concentration strings (R14) to prevent
float-serialization drift.
"""

from __future__ import annotations

import polars as pl

from constants import COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL

PER_CLONOTYPE_SCHEMA: dict[str, pl.DataType] = {
    COL_CLONOTYPE: pl.Utf8,
    "kd": pl.Float64,
    "hillCoefficient": pl.Float64,
    "r2": pl.Float64,
    "affinityClass": pl.Utf8,
    "fitFailureReason": pl.Utf8,
    "kdOutOfRange": pl.Boolean,
}
PER_CLONOTYPE_COLS = list(PER_CLONOTYPE_SCHEMA.keys())


def build_per_clonotype_frame(rows: list[dict]) -> pl.DataFrame:
    """Build the per-clonotype output frame from a list of row dicts.

    Expected keys per row: clonotypeKey, kd (nullable), hillCoefficient (nullable),
    r2 (nullable), affinityClass, fitFailureReason (nullable), kdOutOfRange (nullable).
    """
    if not rows:
        return pl.DataFrame({k: [] for k in PER_CLONOTYPE_SCHEMA}, schema=PER_CLONOTYPE_SCHEMA)
    return pl.DataFrame(rows, schema=PER_CLONOTYPE_SCHEMA)


def flag_kd_out_of_range(
    frame: pl.DataFrame, min_concentration: float, max_concentration: float
) -> pl.DataFrame:
    """R14b: set kdOutOfRange = true when K_D is outside [min_concentration, max_concentration].

    Boundary (kd == min or kd == max) is treated as in-range (closed interval).
    Null K_D rows retain kdOutOfRange = null.
    """
    return frame.with_columns(
        pl.when(pl.col("kd").is_null())
        .then(None)
        .otherwise(
            (pl.col("kd") < min_concentration) | (pl.col("kd") > max_concentration)
        )
        .alias("kdOutOfRange")
    )


def build_mean_bin_frame(signal_frame: pl.DataFrame) -> pl.DataFrame:
    """R14: per-(clonotype, concentration) observed signal (mean_bin or frequency).

    c=0 rows are excluded — they are baseline fixers, not output values.
    Output columns: clonotypeKey, concentrationStr, concentration, meanBin.
    """
    return (
        signal_frame.filter(pl.col(COL_CONC_VAL) != 0)
        .select(
            [
                COL_CLONOTYPE,
                COL_CONC_STR,
                COL_CONC_VAL,
                pl.col("signal").alias("meanBin"),
            ]
        )
    )


FITTED_MEAN_BIN_SCHEMA: dict[str, pl.DataType] = {
    COL_CLONOTYPE: pl.Utf8,
    COL_CONC_STR: pl.Utf8,
    COL_CONC_VAL: pl.Float64,
    "fittedMeanBin": pl.Float64,
}


def build_fitted_mean_bin_frame(fitted_rows: list[dict]) -> pl.DataFrame:
    """R14: per-(clonotype, concentration) fitted signal at experimental concentrations.

    fitted_rows: list of dicts with keys clonotypeKey, concentrationStr, concentration, fittedMeanBin.
    Failed fits contribute no rows (null via absence).
    """
    if not fitted_rows:
        return pl.DataFrame({k: [] for k in FITTED_MEAN_BIN_SCHEMA}, schema=FITTED_MEAN_BIN_SCHEMA)
    return pl.DataFrame(fitted_rows, schema=FITTED_MEAN_BIN_SCHEMA)
