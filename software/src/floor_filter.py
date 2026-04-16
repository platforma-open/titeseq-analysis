"""R8 weight computation + R9 read-floor filtering and insufficient-* classification."""

from __future__ import annotations

import polars as pl

from constants import (
    COL_CLONOTYPE,
    COL_CONC_STR,
    COL_CONC_VAL,
    FitParams,
)
from normalization import CLONOTYPE_READS_AT_CONC, SIGNAL

WEIGHT = "weight"


def apply_floor_and_weights(
    signal_frame: pl.DataFrame, params: FitParams
) -> pl.DataFrame:
    """Apply R8 weights (w_j = clonotype_reads_at_conc) and R9 read floor.

    Rows with `clonotype_reads_at_conc < min_reads_per_concentration` are dropped;
    remaining rows get a `weight` column equal to the per-(clonotype, conc) read sum.
    """
    kept = signal_frame.filter(
        pl.col(CLONOTYPE_READS_AT_CONC) >= params.min_reads_per_concentration
    ).with_columns(pl.col(CLONOTYPE_READS_AT_CONC).cast(pl.Float64).alias(WEIGHT))
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
    counts = non_zero.group_by(COL_CLONOTYPE).agg(pl.len().alias("n_points"))

    all_df = pl.DataFrame({COL_CLONOTYPE: all_clonotypes})
    joined = all_df.join(counts, on=COL_CLONOTYPE, how="left").with_columns(
        pl.col("n_points").fill_null(0)
    )

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
    c0 = filtered.filter(pl.col(COL_CONC_VAL) == 0).select(
        [COL_CLONOTYPE, COL_CONC_STR, SIGNAL, WEIGHT]
    )
    non_c0 = filtered.filter(pl.col(COL_CONC_VAL) != 0)
    return c0, non_c0
