"""R13/R14/R14b: assemble per-clonotype and per-(clonotype, concentration) output frames.

All axis joins are keyed on canonical concentration strings (R14) to prevent
float-serialization drift.
"""

from __future__ import annotations

import polars as pl

from constants import COL_CLONOTYPE, COL_CONC_AM, COL_CONC_STR, COL_CONC_VAL, CONC_AM_SCALE

PER_CLONOTYPE_SCHEMA: dict[str, pl.DataType] = {
    COL_CLONOTYPE: pl.Utf8,
    "kd": pl.Float64,
    "hillCoefficient": pl.Float64,
    "r2": pl.Float64,
    "affinityClass": pl.Utf8,
    "fitFailureReason": pl.Utf8,
    "kdOutOfRange": pl.Boolean,
}


def flag_kd_out_of_range(frame: pl.DataFrame, min_concentration: float, max_concentration: float) -> pl.DataFrame:
    """R14b: set kdOutOfRange = true when K_D is outside [min_concentration, max_concentration].

    Boundary (kd == min or kd == max) is treated as in-range (closed interval).
    Null K_D rows retain kdOutOfRange = null.
    """
    return frame.with_columns(
        pl.when(pl.col("kd").is_null())
        .then(None)
        .otherwise((pl.col("kd") < min_concentration) | (pl.col("kd") > max_concentration))
        .alias("kdOutOfRange")
    )


def add_diagnostic_plot_columns(frame: pl.DataFrame, max_concentration: float) -> pl.DataFrame:
    """R17: append plot-only positions so Failed rows with null K_D render on the scatter.

    kdPlotPosition places null K_D at one decade right of the fitted range on a log axis;
    hillPlotPosition parks null Hill coefficients at 1.0 (centered in the typical [0.5, 2.0] band).
    Both columns are diagnostic — not for reporting — and never null.
    """
    return frame.with_columns(
        pl.when(pl.col("kd").is_null())
        .then(max_concentration * 10.0)
        .otherwise(pl.col("kd"))
        .alias("kdPlotPosition"),
        pl.when(pl.col("hillCoefficient").is_null())
        .then(1.0)
        .otherwise(pl.col("hillCoefficient"))
        .alias("hillPlotPosition"),
    )


def build_mean_bin_frame(signal_frame: pl.DataFrame) -> pl.DataFrame:
    """R14: per-(clonotype, concentration) observed signal (mean_bin or frequency).

    c=0 rows are excluded — they are baseline fixers, not output values.
    Output columns: clonotypeKey, concentrationStr, concentrationAM, concentration, meanBin.
    concentrationAM (attomolar Int64) is the numeric axis key for graph-maker;
    concentrationStr is retained for debugging and backward compatibility.
    """
    return signal_frame.filter(pl.col(COL_CONC_VAL) != 0).select(
        [
            COL_CLONOTYPE,
            COL_CONC_STR,
            (pl.col(COL_CONC_VAL) * CONC_AM_SCALE).round().cast(pl.Int64).alias(COL_CONC_AM),
            COL_CONC_VAL,
            pl.col("signal").alias("meanBin"),
        ]
    )


FITTED_MEAN_BIN_SCHEMA: dict[str, pl.DataType] = {
    COL_CLONOTYPE: pl.Utf8,
    COL_CONC_STR: pl.Utf8,
    COL_CONC_AM: pl.Int64,
    COL_CONC_VAL: pl.Float64,
    "fittedMeanBin": pl.Float64,
}
