"""Pre-fit gating (R6, R8, R9, R9b).

Consolidates everything that runs on the signal frame before Hill fitting:
  - R8 weights + R9 read-floor filtering
  - insufficient_reads / insufficient_points classification
  - c=0 baseline split
  - R6 global baseline B
  - R9b non-monotonic signal (hook effect) detection

Each helper returns a polars frame or scalar; none mutate their inputs.
"""

from __future__ import annotations

import polars as pl

from constants import COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL, FitParams
from normalization import CLONOTYPE_READS_AT_CONC, SIGNAL

WEIGHT = "weight"


def apply_floor_and_weights(signal_frame: pl.DataFrame, params: FitParams) -> pl.DataFrame:
    """Apply R8 weights (w_j = clonotype_reads_at_conc) and R9 read floor.

    Rows with `clonotype_reads_at_conc < min_reads_per_concentration` are dropped;
    remaining rows get a `weight` column equal to the per-(clonotype, conc) read sum.
    """
    kept = signal_frame.filter(pl.col(CLONOTYPE_READS_AT_CONC) >= params.min_reads_per_concentration).with_columns(
        pl.col(CLONOTYPE_READS_AT_CONC).cast(pl.Float64).alias(WEIGHT)
    )
    return kept


def classify_insufficient(
    filtered: pl.DataFrame,
    all_clonotypes: list[str],
    params: FitParams,
) -> pl.DataFrame:
    """R9: mark clonotypes with zero surviving points as insufficient_reads,
    fewer than min_concentration_points as insufficient_points.

    Excludes c=0 points from the count (they are baseline fixers, not fit points).

    Returns DataFrame with columns: clonotypeKey, insufficient_reason (nullable).
    """
    non_zero = filtered.filter(pl.col(COL_CONC_VAL) != 0)
    counts = non_zero.group_by(COL_CLONOTYPE, maintain_order=True).agg(pl.len().alias("n_points"))

    all_df = pl.DataFrame({COL_CLONOTYPE: all_clonotypes})
    joined = all_df.join(counts, on=COL_CLONOTYPE, how="left").with_columns(pl.col("n_points").fill_null(0))

    reason = (
        pl.when(pl.col("n_points") == 0)
        .then(pl.lit("insufficient_reads"))
        .when(pl.col("n_points") < params.min_concentration_points)
        .then(pl.lit("insufficient_points"))
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )
    return joined.with_columns(reason.alias("insufficient_reason")).select(
        [COL_CLONOTYPE, "n_points", "insufficient_reason"]
    )


def split_c0(filtered: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Partition floor-passed points into (c=0 rows for baseline, fit rows)."""
    c0 = filtered.filter(pl.col(COL_CONC_VAL) == 0).select([COL_CLONOTYPE, COL_CONC_STR, SIGNAL, WEIGHT])
    non_c0 = filtered.filter(pl.col(COL_CONC_VAL) != 0)
    return c0, non_c0


def compute_global_baseline(c0_points: pl.DataFrame) -> float | None:
    """R6: arithmetic mean of signal at c=0 across clonotypes that survived R8 floor.

    c0_points is the floor-passed subset of the signal frame where concentration == 0.
    Returns None if the frame is empty (downstream uses 4-param fit).
    """
    if c0_points.height == 0:
        return None
    mean_val = c0_points.select(pl.col(SIGNAL).mean()).item()
    if mean_val is None:
        return None
    return float(mean_val)


def detect_hook_effect(fit_points: pl.DataFrame, bin_mode: bool, params: FitParams) -> pl.DataFrame:
    """R9b: flag a hook only when BOTH spec conditions hold:

      1. (top2_signal - top1_signal) > threshold
      2. (top3_signal - top1_signal) > threshold / 2

    Read-coverage gate per spec requires only top-1 and top-2 reads to clear
    `hook_effect_min_reads` — the top-3 signal comparison is still evaluated, and
    its null (for clonotypes with fewer than three non-zero concentration points)
    collapses the whole expression to False via `fill_null`. Strict `>` is used
    throughout so boundary values do NOT flag.

    Returns DataFrame: clonotypeKey, hook_flag (bool).
    """
    threshold = params.hook_effect_threshold_bin if bin_mode else params.hook_effect_threshold_no_bin
    ranked = fit_points.with_columns(
        pl.col(COL_CONC_VAL).rank(method="ordinal", descending=True).over(COL_CLONOTYPE).alias("rank")
    ).filter(pl.col("rank") <= 3)

    wide = ranked.group_by(COL_CLONOTYPE, maintain_order=True).agg(
        pl.col(SIGNAL).filter(pl.col("rank") == 1).first().alias("top1_signal"),
        pl.col(SIGNAL).filter(pl.col("rank") == 2).first().alias("top2_signal"),
        pl.col(SIGNAL).filter(pl.col("rank") == 3).first().alias("top3_signal"),
        pl.col(CLONOTYPE_READS_AT_CONC).filter(pl.col("rank") == 1).first().alias("top1_reads"),
        pl.col(CLONOTYPE_READS_AT_CONC).filter(pl.col("rank") == 2).first().alias("top2_reads"),
    )

    hook = (
        (pl.col("top1_reads") >= params.hook_effect_min_reads)
        & (pl.col("top2_reads") >= params.hook_effect_min_reads)
        & ((pl.col("top2_signal") - pl.col("top1_signal")) > threshold)
        & ((pl.col("top3_signal") - pl.col("top1_signal")) > threshold / 2)
    )
    return wide.with_columns(hook.fill_null(False).alias("hook_flag")).select([COL_CLONOTYPE, "hook_flag"])
