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

from constants import COL_CLONOTYPE, COL_CONC_LABEL, COL_CONC_STR, COL_CONC_VAL, COL_REPLICATE

PER_CLONOTYPE_SCHEMA: dict[str, pl.DataType] = {
    COL_CLONOTYPE: pl.Utf8,
    "kd": pl.Float64,
    "hillCoefficient": pl.Float64,
    "r2": pl.Float64,
    "affinityClass": pl.Utf8,
    "fitFailureReason": pl.Utf8,
    "kdOutOfRange": pl.Boolean,
    # SE of log10 K_D. Written by "log_fluorescence" only; null under "hill_gate_index".
    "kdLog10Se": pl.Float64,
    # Written by "log_fluorescence" only: "Below range", "In range" or "Above range", null where
    # there is no K_D. Null under "hill_gate_index".
    "kdRangeSide": pl.Utf8,
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


def add_neg_log10_kd(frame: pl.DataFrame) -> pl.DataFrame:
    """Append -log10 of K_D, the form every paper on this assay publishes: larger is tighter.

    The magnitude is the unit the fit solves in and the unit kdLog10Se and kdLog10Sd carry, so the
    value and its error share a scale. Derived, so it is appended after the frame is built rather
    than carried in PER_CLONOTYPE_SCHEMA. Null where there is no K_D, and for a pinned K_D it is
    the pinned value's.
    """
    return frame.with_columns(
        pl.when(pl.col("kd").is_null() | (pl.col("kd") <= 0.0))
        .then(None)
        .otherwise(-pl.col("kd").log10())
        .alias("kdNegLog10")
    )


def add_diagnostic_plot_columns(frame: pl.DataFrame, max_concentration: float) -> pl.DataFrame:
    """R17: append plot-only positions so Failed rows with null Kd render on the scatter.

    kdPlotPosition places null Kd at one decade right of the fitted range on a log axis;
    hillPlotPosition parks null Hill coefficients at -1.0, a non-physical sentinel outside
    the valid (0, ∞) range, so Failed rows pool visually away from well-fitted n≈1
    clonotypes. Both columns are diagnostic — not for reporting — and never null.
    """
    return frame.with_columns(
        pl.when(pl.col("kd").is_null()).then(max_concentration * 10.0).otherwise(pl.col("kd")).alias("kdPlotPosition"),
        pl.when(pl.col("hillCoefficient").is_null())
        .then(-1.0)
        .otherwise(pl.col("hillCoefficient"))
        .alias("hillPlotPosition"),
    )


def build_mean_bin_frame(signal_frame: pl.DataFrame, replicate_keys: list[str] | None = None) -> pl.DataFrame:
    """R14: per-(clonotype, concentration) observed signal (mean_bin or frequency).

    `replicate_keys` names the replicate column when the run has one; the signal then carries it
    as a third key, so the curve of each replicate stays its own.

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
    reps = replicate_keys or []
    return (
        signal_frame.filter(pl.col(COL_CONC_VAL) != 0)
        .select(
            [
                COL_CLONOTYPE,
                *reps,
                COL_CONC_STR,
                pl.col("signal").alias("meanBin"),
            ]
        )
        .sort([COL_CLONOTYPE, *reps, COL_CONC_STR])
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


CONCENTRATION_LABEL_SCHEMA: dict[str, pl.DataType] = {
    COL_CONC_STR: pl.Utf8,
    COL_CONC_LABEL: pl.Utf8,
}


def build_concentration_label_frame(concentration_value: pl.DataFrame) -> pl.DataFrame:
    """Each non-zero concentration as a category label that sorts as text in numeric order.

    The rank, zero-padded, then the concentration: "01 · 3.16e-10". The category graphs group by
    it, since they sort a String axis as text and ignore a stored order once they facet. Two strings
    of one value rank in text order.
    """
    width = len(str(concentration_value.height))
    rank = pl.int_range(1, pl.len() + 1).cast(pl.Utf8).str.zfill(width)
    return (
        concentration_value.sort([COL_CONC_VAL, COL_CONC_STR])
        .select(COL_CONC_STR, (rank + pl.lit(" · ") + pl.col(COL_CONC_STR)).alias(COL_CONC_LABEL))
        .sort(COL_CONC_STR)
    )


COL_SERIES = "series"
SERIES_OBSERVED = "Observed"
SERIES_FITTED = "Fitted"


def _curved_signal(mean_bin: pl.DataFrame, fitted: pl.DataFrame) -> tuple[pl.DataFrame, list[str]]:
    """The observed signal of the fits that drew a curve, and the key columns the two frames share.

    Only fits with a curve count, so a mean over replicates reads the same replicates for the points
    and for the curve.
    """
    keys = [c for c in mean_bin.columns if c != "meanBin" and c in fitted.columns]
    curved = fitted.select([c for c in keys if c != COL_CONC_STR]).unique()
    return mean_bin.join(curved, on=curved.columns, how="semi"), keys


def build_signal_series_frame(mean_bin: pl.DataFrame, fitted: pl.DataFrame) -> pl.DataFrame:
    """The observed signal and the fitted curve stacked in one value column, told apart by series.

    A category graph takes one value column: stacked, the two draw as two lines in each panel.
    Keeps every key column the two frames share (clonotype, replicate where present, concentration).
    """
    observed, keys = _curved_signal(mean_bin, fitted)
    points = observed.select(*keys, pl.lit(SERIES_OBSERVED).alias(COL_SERIES), pl.col("meanBin").alias("value"))
    curve = fitted.select(*keys, pl.lit(SERIES_FITTED).alias(COL_SERIES), pl.col("fittedMeanBin").alias("value"))
    return pl.concat([points, curve]).sort([*keys, COL_SERIES])


def build_combined_signal_frame(mean_bin: pl.DataFrame, fitted: pl.DataFrame) -> pl.DataFrame:
    """Per variant and concentration: the observed mean over replicates, its standard error, and the mean curve.

    The standard error is the sample SD over sqrt(n), null below two replicates. Every concentration
    of the curve is kept; the observed columns are null where nothing was measured. The observed mean
    reads the replicates measured at that concentration; the curve mean reads every replicate with a
    curve, since the pipeline writes each curve on the union of the replicates' concentrations. The
    column names differ from the per-replicate frames': the workflow keys each imported column by its
    name.
    """
    observed, keys = _curved_signal(mean_bin, fitted)
    by = [c for c in keys if c != COL_REPLICATE]
    points = observed.group_by(by).agg(
        pl.col("meanBin").mean().alias("combinedMeanBin"),
        (pl.col("meanBin").std(ddof=1) / pl.col("meanBin").count().sqrt()).alias("combinedMeanBinSe"),
    )
    return (
        fitted.group_by(by)
        .agg(pl.col("fittedMeanBin").mean().alias("combinedFittedMeanBin"))
        .join(points, on=by, how="full", coalesce=True)
        .select(*by, "combinedMeanBin", "combinedMeanBinSe", "combinedFittedMeanBin")
        .sort(by)
    )
