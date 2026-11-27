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

from constants import COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL, COL_REPLICATE

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


def add_log10_kd(frame: pl.DataFrame) -> pl.DataFrame:
    """Append log10 of K_D and -log10 of K_D, the form every paper on this assay publishes: larger
    is tighter.

    The magnitude is the unit the fit solves in and the unit kdLog10Se and kdLog10Sd carry, so the
    value and its error share a scale. Derived, so they are appended after the frame is built rather
    than carried in PER_CLONOTYPE_SCHEMA. Null where there is no K_D, and for a pinned K_D they are
    the pinned value's.
    """
    log10_kd = pl.when(pl.col("kd").is_null() | (pl.col("kd") <= 0.0)).then(None).otherwise(pl.col("kd").log10())
    return frame.with_columns(log10_kd.alias("kdLog10"), (-log10_kd).alias("kdNegLog10"))


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


COL_AA_KEY = "aaKey"


def curve_variants(nt_to_aa: pl.DataFrame, fitted: pl.DataFrame) -> pl.DataFrame:
    """Synonymous mode: [clonotypeKey, aaKey] for the nt variants that drew a fitted curve.

    `nt_to_aa` is [clonotypeKey, aaKey], one amino-acid variant per nt variant.
    """
    curved = fitted.select(pl.col(COL_CLONOTYPE).cast(pl.Utf8)).unique()
    link = nt_to_aa.select(pl.col(COL_CLONOTYPE).cast(pl.Utf8), pl.col(COL_AA_KEY).cast(pl.Utf8))
    return link.join(curved, on=COL_CLONOTYPE, how="semi").sort(COL_CLONOTYPE)


# A standard error of log10 Kd,app above this many decades gets no range: 10^SE is not a number.
_MARK_SE_MAX = 100.0
# The Kd,app ± SE band is one line that crosses it this many times per decade of x, and never fewer
# than the second number, see `_kd_mark_rows`. A chart shows a decade about 100 px wide and can be
# enlarged twofold: the line, 6 px wide, must still cover the gap between two crossings.
_BAND_STROKES_PER_DECADE = 64
_BAND_STROKES_MIN = 8

_CURVE_VALUES = [
    "curveFittedMeanBin",
    "signalConcentration",
    "curveCombinedFittedMeanBin",
    "curveExtent",
    "kdPoint",
    "kdReplicatePoint",
    "kdSeRange",
]


def _height_at(points: pl.DataFrame, at: pl.DataFrame, by: list[str]) -> pl.DataFrame:
    """`at` ([*by, kd]) with `height`: each curve of `points` ([*by, x, y]) read at its `kd`.

    Straight between the two points either side, on a log x. Beyond the curve's ends, the height of
    the nearest end. A group of `at` without a curve is dropped.
    """
    ordered = points.sort([*by, "x"])
    segments = ordered.with_columns(
        pl.col("x").shift(-1).over(by).fill_null(pl.col("x")).alias("x1"),
        pl.col("y").shift(-1).over(by).fill_null(pl.col("y")).alias("y1"),
        pl.col("x").min().over(by).alias("x_low"),
        pl.col("x").max().over(by).alias("x_high"),
    )
    x, x1, y, y1 = pl.col("x"), pl.col("x1"), pl.col("y"), pl.col("y1")
    inside = pl.col("kd").clip(pl.col("x_low"), pl.col("x_high"))
    share = pl.when(x1 > x).then((inside.log10() - x.log10()) / (x1.log10() - x.log10())).otherwise(0.0)
    return (
        segments.join(at, on=by, how="inner")
        .filter((x <= inside) & (inside <= x1))
        .with_columns((y + (y1 - y) * share).alias("height"))
        .unique(subset=by, keep="first", maintain_order=True)
        .select(*by, "kd", "height")
    )


def _kd_mark_rows(
    frame: pl.DataFrame,
    kd: pl.DataFrame,
    replicate_kd: pl.DataFrame | None,
    key: str,
    x_range: tuple[float, float] | None,
) -> pl.DataFrame:
    """The Kd,app marks of each variant, as points of three more curves.

    `kdPoint`: the variant's Kd,app, on its mean curve. `kdSeRange`: a band from Kd,app / 10^SE to
    Kd,app × 10^SE, over the height of the variant's curves. A chart draws lines only, so the band
    is one line that runs between the bottom and the top `_BAND_STROKES_PER_DECADE` times a decade,
    at even steps of log x; drawn wider than a step, it reads as filled. The band stops at the ends of `x_range`,
    the lowest and highest non-zero concentrations, or at Kd,app where that lies outside them: it
    never widens the chart. Both sit on the variant's first replicate alone: a chart draws a curve
    once per replicate, and a band drawn several times darkens. `kdReplicatePoint`: each replicate's own
    Kd,app from `replicate_kd` ([key, replicate, kd]), on that replicate's curve. The points have
    concentration keys of their own, and `curveExtent` holds each point's height. A variant without
    Kd,app has no `kdPoint`; one without SE, or with one beyond `_MARK_SE_MAX`, no band.
    """
    schema = {
        key: pl.Utf8,
        COL_REPLICATE: pl.Utf8,
        COL_CONC_STR: pl.Utf8,
        "signalConcentration": pl.Float64,
        "curveExtent": pl.Float64,
        "kdPoint": pl.Float64,
        "kdReplicatePoint": pl.Float64,
        "kdSeRange": pl.Float64,
    }
    positive = pl.col("kd").is_not_null() & (pl.col("kd") > 0)
    mean_curve = frame.select(
        key, pl.col("signalConcentration").alias("x"), pl.col("curveCombinedFittedMeanBin").alias("y")
    ).unique()
    spans = frame.group_by(key).agg(
        pl.col(COL_REPLICATE).min(),
        pl.col("curveFittedMeanBin").min().alias("low"),
        pl.col("curveFittedMeanBin").max().alias("high"),
    )
    kd_value, se = pl.col("kd"), pl.col("kdLog10Se")
    x_min, x_max = (kd_value, kd_value) if x_range is None else (pl.lit(x_range[0]), pl.lit(x_range[1]))
    combined = (
        _height_at(mean_curve, kd.select(pl.col(key).cast(pl.Utf8), "kd").filter(positive), [key])
        .join(kd.select(pl.col(key).cast(pl.Utf8), "kdLog10Se"), on=key, how="left")
        .join(spans, on=key, how="left")
        .with_columns(
            pl.max_horizontal(kd_value / 10.0**se, pl.min_horizontal(x_min, kd_value)).alias("lo"),
            pl.min_horizontal(kd_value * 10.0**se, pl.max_horizontal(x_max, kd_value)).alias("hi"),
            (se.is_not_null() & (se > 0) & (se <= _MARK_SE_MAX)).alias("ranged"),
        )
    )
    stroke, strokes = pl.col("stroke"), pl.col("strokes")
    decades = (pl.col("hi") / pl.col("lo")).log10()
    band = (
        combined.filter("ranged")
        .with_columns(
            (decades * _BAND_STROKES_PER_DECADE)
            .ceil()
            .clip(lower_bound=_BAND_STROKES_MIN)
            .cast(pl.Int64)
            .alias("strokes")
        )
        .with_columns(pl.int_ranges(0, strokes + 1).alias("stroke"))
        .explode("stroke")
        .with_columns(
            (pl.col("lo") * (pl.col("hi") / pl.col("lo")) ** (stroke / strokes)).alias("x"),
            pl.when(stroke % 2 == 0).then(pl.col("low")).otherwise(pl.col("high")).alias("height"),
            pl.format("kd-se-{}", stroke.cast(pl.Utf8).str.zfill(4)).alias("name"),
        )
    )
    # (column, the rows that have the point, its concentration key, x)
    points = [
        ("kdPoint", combined, pl.lit("kd-point"), kd_value),
        ("kdSeRange", band, pl.col("name"), pl.col("x")),
    ]
    if replicate_kd is not None:
        by = [key, COL_REPLICATE]
        replicate_curves = frame.select(
            *by, pl.col("signalConcentration").alias("x"), pl.col("curveFittedMeanBin").alias("y")
        )
        own = replicate_kd.select(pl.col(key).cast(pl.Utf8), pl.col(COL_REPLICATE).cast(pl.Utf8), "kd")
        heights = _height_at(replicate_curves, own.filter(positive), by)
        points.append(("kdReplicatePoint", heights, pl.lit("kd-replicate"), kd_value))
    frames = [
        rows.select(
            key,
            COL_REPLICATE,
            name.alias(COL_CONC_STR),
            x.alias("signalConcentration"),
            pl.col("height").alias("curveExtent"),
            pl.col("height").alias(column),
        )
        for column, rows, name, x in points
    ]
    return pl.concat([pl.DataFrame(schema=schema), *frames], how="diagonal").select(list(schema))


def _curve_frame(
    curves: pl.DataFrame,
    concentration_value: pl.DataFrame,
    key: str,
    kd: pl.DataFrame | None,
    replicate_kd: pl.DataFrame | None,
) -> pl.DataFrame:
    """One chart's worth of columns per variant, all on [key, replicate, concentrationStr].

    `curves` is [key, replicate, concentrationStr, curveFittedMeanBin]. Adds the mean curve over
    replicates (repeated on every replicate's row), the numeric concentration for a scatterplot's
    x, and from `kd` ([key, kd, kdLog10Se, ...]) and `replicate_kd` each variant's Kd,app marks,
    see `_kd_mark_rows`. `curveExtent` is set on every row: the mean curve on a curve's row, the
    mark's height on a mark's.
    """
    by = [key, COL_REPLICATE, COL_CONC_STR]
    combined = curves.group_by([key, COL_CONC_STR]).agg(
        pl.col("curveFittedMeanBin").mean().alias("curveCombinedFittedMeanBin")
    )
    measured = concentration_value.filter(pl.col(COL_CONC_VAL) > 0)[COL_CONC_VAL]
    x_range = None if measured.is_empty() else (measured.min(), measured.max())
    frame = (
        curves.join(combined, on=[key, COL_CONC_STR], how="left")
        .join(
            concentration_value.select(COL_CONC_STR, pl.col(COL_CONC_VAL).alias("signalConcentration")),
            on=COL_CONC_STR,
            how="inner",
        )
        .with_columns(pl.col("curveCombinedFittedMeanBin").alias("curveExtent"))
    )
    marks = _kd_mark_rows(
        frame.filter(pl.col("signalConcentration") > 0),
        pl.DataFrame(schema={key: pl.Utf8, "kd": pl.Float64, "kdLog10Se": pl.Float64}) if kd is None else kd,
        replicate_kd,
        key,
        x_range,
    )
    return pl.concat([frame, marks], how="diagonal").select(*by, *_CURVE_VALUES).sort(by, nulls_last=False)


def _with_replicate(frame: pl.DataFrame) -> pl.DataFrame:
    """`frame` with a text replicate column: null on a run without replicates."""
    if COL_REPLICATE not in frame.columns:
        return frame.with_columns(pl.lit(None, dtype=pl.Utf8).alias(COL_REPLICATE))
    return frame.with_columns(pl.col(COL_REPLICATE).cast(pl.Utf8))


def build_variant_curve_frame(
    fitted: pl.DataFrame,
    concentration_value: pl.DataFrame,
    per_clonotype: pl.DataFrame | None = None,
    per_replicate: pl.DataFrame | None = None,
) -> pl.DataFrame:
    """The variant panel's chart of each fitted variant, beside its fitted curve: see `_curve_frame`.

    [clonotypeKey, replicate, concentrationStr] with the mean curve over replicates, the numeric
    concentration and the Kd,app marks of `per_clonotype` and, with replicates, of `per_replicate`
    ([clonotypeKey, replicate, kd, ...]). The fitted curve itself is `fitted`, on the same axes, and
    is not repeated. A run without replicates has no replicate column.
    """
    replicated = COL_REPLICATE in fitted.columns
    curves = _with_replicate(fitted).select(
        pl.col(COL_CLONOTYPE).cast(pl.Utf8),
        COL_REPLICATE,
        COL_CONC_STR,
        pl.col("fittedMeanBin").alias("curveFittedMeanBin"),
    )
    out = _curve_frame(
        curves, concentration_value, COL_CLONOTYPE, per_clonotype, per_replicate if replicated else None
    ).drop("curveFittedMeanBin")
    return out if replicated else out.drop(COL_REPLICATE)


_AA_CURVE_NAMES = {
    "curveFittedMeanBin": "aaFittedMeanBin",
    "signalConcentration": "aaSignalConcentration",
    "curveCombinedFittedMeanBin": "aaCombinedFittedMeanBin",
    "curveExtent": "aaCurveExtent",
    "kdPoint": "aaKdPoint",
    "kdReplicatePoint": "aaKdReplicatePoint",
    "kdSeRange": "aaKdSeRange",
}


def build_aa_curve_frame(
    fitted: pl.DataFrame,
    nt_to_aa: pl.DataFrame,
    concentration_value: pl.DataFrame,
    per_variant: pl.DataFrame | None = None,
    per_replicate: pl.DataFrame | None = None,
) -> pl.DataFrame:
    """Synonymous mode with replicates: each amino-acid variant's curves, made from its nt variants' fits.

    Nothing is fitted at the amino-acid grain. Per (amino-acid variant, replicate, concentration): the
    mean over its nt variants of the fitted curve. `nt_to_aa` is [clonotypeKey, aaKey], one
    amino-acid variant per nt variant. The rest is `_curve_frame`'s, under the amino-acid column
    names; `per_variant` ([aaKey, kd, kdLog10Se, ...]) and `per_replicate` ([aaKey, replicate, kd,
    ...]) give the Kd,app marks.
    """
    by = [COL_AA_KEY, COL_REPLICATE, COL_CONC_STR]
    # Only a run where nothing was fitted lacks the replicate column.
    fitted = _with_replicate(fitted)
    link = nt_to_aa.select(pl.col(COL_CLONOTYPE).cast(pl.Utf8), pl.col(COL_AA_KEY).cast(pl.Utf8))
    curved = fitted.with_columns(pl.col(COL_CLONOTYPE).cast(pl.Utf8)).join(link, on=COL_CLONOTYPE, how="inner")
    curves = curved.group_by(by).agg(pl.col("fittedMeanBin").mean().alias("curveFittedMeanBin"))
    return _curve_frame(curves, concentration_value, COL_AA_KEY, per_variant, per_replicate).rename(_AA_CURVE_NAMES)
