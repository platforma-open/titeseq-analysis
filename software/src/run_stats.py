"""Run-level numbers for the summary pages: classes, walls, concentrations, parents, read depth."""

from __future__ import annotations

import math
import statistics

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


def _wall_side(kd: float | None, out_of_range: bool | None, lo: float | None, hi: float | None) -> str | None:
    if not out_of_range or kd is None or kd <= 0 or lo is None or hi is None:
        return None
    return "below" if math.log10(kd) <= (math.log10(lo) + math.log10(hi)) / 2 else "above"


def run_stats(
    frame: pl.DataFrame,
    *,
    grain: str,
    conc_bounds: tuple[float | None, float | None],
    concentrations: int,
    replicates: list[str],
    parents: list[dict] | None = None,
) -> dict:
    """The JSON-ready run statistics of the reported frame: one row per reported variant, with kd,
    affinityClass, kdOutOfRange and fitReads."""
    lo, hi = conc_bounds
    classes: dict[str, int] = {}
    walls = {"below": 0, "above": 0}
    reads_by_class: dict[str, list[int]] = {}
    for row in frame.select("kd", "affinityClass", "kdOutOfRange", COL_FIT_READS).iter_rows(named=True):
        cls = row["affinityClass"] or "Failed"
        classes[cls] = classes.get(cls, 0) + 1
        side = _wall_side(row["kd"], row["kdOutOfRange"], lo, hi)
        if side is not None:
            walls[side] += 1
        reads_by_class.setdefault(cls, []).append(row[COL_FIT_READS])
    stats: dict = {
        "grain": grain,
        "variants": frame.height,
        "classes": classes,
        "walls": walls,
        "concentrations": {
            "count": concentrations,
            "log10Lowest": math.log10(lo) if lo else None,
            "log10Highest": math.log10(hi) if hi else None,
        },
        "replicates": replicates,
        "medianReadsByClass": {cls: statistics.median(v) for cls, v in sorted(reads_by_class.items())},
    }
    if parents is not None:
        stats["parents"] = [
            {
                **p,
                "tightRoom": None if p["log10Kd"] is None or not lo else p["log10Kd"] - math.log10(lo),
                "weakRoom": None if p["log10Kd"] is None or not hi else math.log10(hi) - p["log10Kd"],
            }
            for p in parents
        ]
    return stats
