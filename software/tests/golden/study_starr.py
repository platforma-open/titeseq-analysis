"""Starr 2020 and 2022 as the block sees them: the authors' barcode counts per sample, one library at a
time, four gates with their sort fractions, run through `pipeline.run` under the mean bin fit.

Each barcode is a variant. A fixed sample of the barcodes is fitted (`paper_starr.SAMPLE_ENDINGS`); the
others are summed into one filler variant per sample, so gate depths and mean bins are the full
library's.
"""

from __future__ import annotations

from functools import cache

import polars as pl
from golden_support import cached_run, study_file
from paper_starr import FILLER, LIBRARIES, SAMPLE_ENDINGS

from constants import FitParams

# The block's defaults for the mean bin fit (model/src/fitSettings.ts HILL_GATE_INDEX_DEFAULTS).
PARAMS = FitParams(equation="hill_gate_index")


def _file(study: str, path: str) -> str:
    return str(study_file(study, path))


def _sampled(barcode: pl.Expr) -> pl.Expr:
    # By the barcode's own last bases, not a hash: polars does not keep its hash stable across releases.
    return barcode.str.slice(-2).is_in(SAMPLE_ENDINGS)


@cache
def reads(case: str) -> pl.DataFrame:
    lib = LIBRARIES[case]
    meta = pl.read_csv(_file(lib.study, "upload/metadata.csv")).filter(pl.col("Library") == lib.library)
    counts = pl.read_csv(_file(lib.study, "upload/counts-long.csv")).join(
        meta.select("Sample"), on="Sample", how="semi"
    )
    kept = counts.filter(_sampled(pl.col("Barcode")))
    filler = (
        counts.filter(~_sampled(pl.col("Barcode")))
        .group_by("Sample")
        .agg(pl.col("Count").sum())
        .select("Sample", pl.lit(FILLER).alias("Barcode"), "Count")
    )
    return (
        pl.concat([kept, filler])
        .join(meta, on="Sample")
        .select(
            pl.col("Barcode").alias("clonotypeKey"),
            pl.col("Sample").alias("sampleId"),
            pl.col("ACE2 concentration (M)").cast(pl.Utf8).alias("concentrationStr"),
            pl.col("ACE2 concentration (M)").cast(pl.Float64).alias("concentration"),
            pl.col("FACS bin").cast(pl.Int64).alias("bin"),
            pl.col("Count").cast(pl.Int64).alias("reads"),
            pl.col("Sort fraction").alias("sort_fraction"),
        )
        .sort("sampleId", "clonotypeKey")
    )


def run_block(case: str, params: FitParams = PARAMS) -> dict[str, pl.DataFrame]:
    from pipeline import run

    lib = LIBRARIES[case]
    return cached_run(
        case,
        lambda: run(reads(case), params=params, sort_fraction_col="sort_fraction"),
        params,
        inputs=[_file(lib.study, "upload/metadata.csv"), _file(lib.study, "upload/counts-long.csv")],
    )


@cache
def published(case: str) -> dict[str, float]:
    """barcode -> the authors' log10 Ka (= -log10 Kd), for the sampled barcodes they reported one for."""
    lib = LIBRARIES[case]
    key = pl.read_csv(_file(lib.study, lib.answer_key), null_values="NA", infer_schema_length=10000)
    key = key.filter((pl.col("library") == lib.key_library) & pl.col("log10Ka").is_not_null())
    key = key.filter(_sampled(pl.col("barcode")))
    return dict(key.select("barcode", "log10Ka").iter_rows())


def clip(case: str, neg_log10_kd: float) -> float:
    lib = LIBRARIES[case]
    return min(max(neg_log10_kd, lib.neg_log10_kd_lower), lib.neg_log10_kd_upper)


def in_range(case: str, neg_log10_kd: float) -> bool:
    lib = LIBRARIES[case]
    return lib.neg_log10_kd_lower + 0.01 < neg_log10_kd < lib.neg_log10_kd_upper - 0.01
