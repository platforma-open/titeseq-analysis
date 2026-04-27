"""R13/R14/R14b: assemble per-clonotype and per-(clonotype, concentration) output frames.

Output TSVs carry only the canonical `concentrationStr` column for the
concentration axis (R14: parse-once invariant). The Tengo workflow wraps that
column directly as a String axis on the output PColumns. See
`docs/investigations/concentration-axis-spec-realignment.md` for the spec
deviation rationale (spec calls for Float axis; SDK gates axis types to
Int|Long|String — Graph Maker renders the String axis categorically, ordered
by parsed numeric value).
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


def flag_kd_out_of_range(frame: pl.DataFrame, min_concentration: float, max_concentration: float) -> pl.DataFrame:
    """R14b: set kdOutOfRange = true when Kd is outside [min_concentration, max_concentration].

    Boundary (kd == min or kd == max) is treated as in-range (closed interval).
    Null Kd rows retain kdOutOfRange = null.
    """
    return frame.with_columns(
        pl.when(pl.col("kd").is_null())
        .then(None)
        .otherwise((pl.col("kd") < min_concentration) | (pl.col("kd") > max_concentration))
        .alias("kdOutOfRange")
    )


def add_diagnostic_plot_columns(frame: pl.DataFrame, max_concentration: float) -> pl.DataFrame:
    """R17: append plot-only positions so Failed rows with null Kd render on the scatter.

    kdPlotPosition places null Kd at one decade right of the fitted range on a log axis;
    hillPlotPosition parks null Hill coefficients at -1.0, a non-physical sentinel outside
    the valid (0, ∞) range, so Failed rows pool visually away from well-fitted n≈1
    clonotypes. Both columns are diagnostic — not for reporting — and never null.
    """
    return frame.with_columns(
        pl.when(pl.col("kd").is_null())
        .then(max_concentration * 10.0)
        .otherwise(pl.col("kd"))
        .alias("kdPlotPosition"),
        pl.when(pl.col("hillCoefficient").is_null())
        .then(-1.0)
        .otherwise(pl.col("hillCoefficient"))
        .alias("hillPlotPosition"),
    )


def build_mean_bin_frame(signal_frame: pl.DataFrame) -> pl.DataFrame:
    """R14: per-(clonotype, concentration) observed signal (mean_bin or frequency).

    c=0 rows are excluded — they are baseline fixers, not output values, and on a
    log-scale axis log(0) = −∞ would break Graph Maker rendering anyway.
    Output columns: clonotypeKey, concentrationStr, meanBin. The workflow wraps
    concentrationStr as the canonical String axis directly.
    """
    return signal_frame.filter(pl.col(COL_CONC_VAL) != 0).select(
        [
            COL_CLONOTYPE,
            COL_CONC_STR,
            pl.col("signal").alias("meanBin"),
        ]
    )


FITTED_MEAN_BIN_SCHEMA: dict[str, pl.DataType] = {
    COL_CLONOTYPE: pl.Utf8,
    COL_CONC_STR: pl.Utf8,
    "fittedMeanBin": pl.Float64,
}
