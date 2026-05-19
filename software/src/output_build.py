"""R13/R14/R14b: assemble per-clonotype and per-(clonotype, concentration) output frames.

Output TSVs carry the canonical `concentrationStr` column on the concentration
axis (R14 parse-once invariant). The Tengo workflow wraps it as a String axis
on the output PColumns. Log-scale graph rendering reads numeric values from a
separate `concentrationValue` PColumn (axes `[concentration:String]`, valueType
Double). Spec calls for a Float axis but the SDK regex at
`core/platforma/sdk/workflow-tengo/src/pt/util.lib.tengo:352` limits axis types
to `Int|Long|String`.
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

    Excludes c=0 rows — they fix the baseline, not output values, and log(0) =
    −∞ would break log-scale rendering. Output columns: clonotypeKey,
    concentrationStr, meanBin. The numeric source for log-scale graphs lives in
    build_concentration_value_frame.

    Sorted by (clonotypeKey, concentrationStr) so the TSV is byte-stable across
    runs. Upstream `signal_frame` comes from `polars.group_by(...).agg(...)` in
    `normalize()`, which does not guarantee row order — leaving this unsorted
    let identical inputs produce TSVs whose row order varied run-to-run, which
    propagated through xsv.importFile to a different Parquet content hash and
    triggered CIDConflictError on re-derivation.
    """
    return (
        signal_frame.filter(pl.col(COL_CONC_VAL) != 0)
        .select(
            [
                COL_CLONOTYPE,
                COL_CONC_STR,
                pl.col("signal").alias("meanBin"),
            ]
        )
        .sort([COL_CLONOTYPE, COL_CONC_STR])
    )


def build_concentration_value_frame(signal_frame: pl.DataFrame) -> pl.DataFrame:
    """Sidecar PColumn for the Titration Curves X-axis.

    Axes `[concentration:String]`, valueType Double. One row per unique non-zero
    concentration. Graph Maker reads the Double for log-scale X positions while
    joining to meanBin / fittedMeanBin on the shared String axis. R14: each
    concentrationStr maps to the Float64 already parsed by
    canonicalize_concentration — never re-parsed.
    """
    return (
        signal_frame.filter(pl.col(COL_CONC_VAL) != 0)
        .select([COL_CONC_STR, COL_CONC_VAL])
        .unique()
        .sort([COL_CONC_STR, COL_CONC_VAL])
    )


FITTED_MEAN_BIN_SCHEMA: dict[str, pl.DataType] = {
    COL_CLONOTYPE: pl.Utf8,
    COL_CONC_STR: pl.Utf8,
    "fittedMeanBin": pl.Float64,
}


CONCENTRATION_VALUE_SCHEMA: dict[str, pl.DataType] = {
    COL_CONC_STR: pl.Utf8,
    COL_CONC_VAL: pl.Float64,
}
