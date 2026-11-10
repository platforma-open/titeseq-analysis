"""The reads each variant's fit used, for the read depth columns."""

from __future__ import annotations

import polars as pl

COL_FIT_READS = "fitReads"


def fit_reads(reads: pl.DataFrame, key: str, reads_col: str) -> pl.DataFrame:
    """[key, fitReads]: each key's total of `reads_col` over the rows of `reads`."""
    return (
        reads.group_by(pl.col(key).cast(pl.Utf8))
        .agg(pl.col(reads_col).fill_null(0).sum().cast(pl.Int64).alias(COL_FIT_READS))
        .sort(key)
    )


def with_fit_reads(frame: pl.DataFrame, key: str, totals: dict[str, int]) -> pl.DataFrame:
    """`frame` with a fitReads column read from `totals` by `key`; 0 where a key has none."""
    keys = frame[key].cast(pl.Utf8).to_list()
    return frame.with_columns(pl.Series(COL_FIT_READS, [totals.get(k, 0) for k in keys], dtype=pl.Int64))
