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
    )
    top2 = ranked.filter(pl.col("rank") <= 2)

    pivot_parts = []
    for rank_value, label in ((1, "top1"), (2, "top2")):
        part = (
            top2.filter(pl.col("rank") == rank_value)
            .select([COL_CLONOTYPE, SIGNAL, CLONOTYPE_READS_AT_CONC])
            .rename(
                {SIGNAL: f"{label}_signal", CLONOTYPE_READS_AT_CONC: f"{label}_reads"}
            )
        )
        pivot_parts.append(part)
    wide = pivot_parts[0].join(pivot_parts[1], on=COL_CLONOTYPE, how="full", coalesce=True)

    hook = (
        (pl.col("top1_reads") >= params.hook_effect_min_reads)
        & (pl.col("top2_reads") >= params.hook_effect_min_reads)
        & ((pl.col("top2_signal") - pl.col("top1_signal")) > threshold)
    )
    return wide.with_columns(hook.fill_null(False).alias("hook_flag")).select(
        [COL_CLONOTYPE, "hook_flag"]
    )
