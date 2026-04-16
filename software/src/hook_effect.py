"""R9b: pre-fit non-monotonic signal detection at the top two concentrations."""

from __future__ import annotations

import polars as pl

from constants import COL_CLONOTYPE, COL_CONC_VAL, FitParams
from normalization import CLONOTYPE_READS_AT_CONC, SIGNAL


def detect_hook_effect(
    fit_points: pl.DataFrame, bin_mode: bool, params: FitParams
) -> pl.DataFrame:
    """R9b: signal drop from conc rank-2 to conc rank-1 that exceeds threshold flags a hook.

    Both top-2 and top-1 concentration points must have >= hook_effect_min_reads.
    Uses strict `>` for the drop comparison (drop == threshold does NOT flag).

    Returns DataFrame: clonotypeKey, hook_flag (bool).
    """
    threshold = (
        params.hook_effect_threshold_bin if bin_mode else params.hook_effect_threshold_no_bin
    )
    ranked = fit_points.with_columns(
        pl.col(COL_CONC_VAL)
        .rank(method="ordinal", descending=True)
        .over(COL_CLONOTYPE)
        .alias("rank")
    ).filter(pl.col("rank") <= 2)

    wide = ranked.group_by(COL_CLONOTYPE).agg(
        pl.col(SIGNAL).filter(pl.col("rank") == 1).first().alias("top1_signal"),
        pl.col(SIGNAL).filter(pl.col("rank") == 2).first().alias("top2_signal"),
        pl.col(CLONOTYPE_READS_AT_CONC).filter(pl.col("rank") == 1).first().alias("top1_reads"),
        pl.col(CLONOTYPE_READS_AT_CONC).filter(pl.col("rank") == 2).first().alias("top2_reads"),
    )

    hook = (
        (pl.col("top1_reads") >= params.hook_effect_min_reads)
        & (pl.col("top2_reads") >= params.hook_effect_min_reads)
        & ((pl.col("top2_signal") - pl.col("top1_signal")) > threshold)
    )
    return wide.with_columns(hook.fill_null(False).alias("hook_flag")).select(
        [COL_CLONOTYPE, "hook_flag"]
    )
