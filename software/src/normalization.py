"""Signal computation (R7 mean bin, R7b frequency) — polars-native, no per-row Python."""

from __future__ import annotations

import polars as pl

from constants import COL_BIN, COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL, COL_READS

MEAN_BIN = "mean_bin"
SIGNAL = "signal"  # unified name: mean_bin in bin mode, frequency in no-bin mode
FREQ = "freq"
TOTAL_READS_AT_CONC = "total_reads_at_conc"
CLONOTYPE_READS_AT_CONC = "clonotype_reads_at_conc"


def compute_mean_bin(reads: pl.DataFrame) -> pl.DataFrame:
    """R7: mean_bin_c = Σ_b (b · freq_cb) / Σ_b freq_cb, where freq_cb = reads_cb / depth_b.

    Groups by (clonotype, concentration). Per-bin frequency uses per-sample depth
    (total reads in that bin at that concentration). A bin-at-concentration with zero
    depth contributes nothing to either numerator or denominator.

    Returns long frame: clonotypeKey, concentrationStr, concentration, mean_bin, clonotype_reads_at_conc.
    """
    with_freq = reads.with_columns(
        pl.col(COL_READS).sum().over([COL_CONC_STR, COL_BIN]).alias("depth"),
    ).with_columns(pl.when(pl.col("depth") > 0).then(pl.col(COL_READS) / pl.col("depth")).otherwise(0.0).alias(FREQ))
    return (
        with_freq.group_by([COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL])
        .agg(
            (pl.col(COL_BIN).cast(pl.Float64) * pl.col(FREQ)).sum().alias("num"),
            pl.col(FREQ).sum().alias("den"),
            pl.col(COL_READS).sum().alias(CLONOTYPE_READS_AT_CONC),
        )
        .with_columns(pl.when(pl.col("den") > 0).then(pl.col("num") / pl.col("den")).otherwise(None).alias(MEAN_BIN))
        .drop(["num", "den"])
    )


def compute_frequency_signal(reads: pl.DataFrame) -> pl.DataFrame:
    """R7b: no-bin mode signal = reads_clonotype_at_conc / total_reads_at_conc.

    Returns long frame: clonotypeKey, concentrationStr, concentration, signal, clonotype_reads_at_conc.
    """
    per_clonotype = reads.group_by([COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL]).agg(
        pl.col(COL_READS).sum().alias(CLONOTYPE_READS_AT_CONC)
    )
    return (
        per_clonotype.with_columns(
            pl.col(CLONOTYPE_READS_AT_CONC).sum().over(COL_CONC_STR).alias(TOTAL_READS_AT_CONC),
        )
        .with_columns(
            pl.when(pl.col(TOTAL_READS_AT_CONC) > 0)
            .then(pl.col(CLONOTYPE_READS_AT_CONC) / pl.col(TOTAL_READS_AT_CONC))
            .otherwise(None)
            .alias(SIGNAL)
        )
        .drop(TOTAL_READS_AT_CONC)
    )


def normalize(reads: pl.DataFrame, bin_mode: bool) -> pl.DataFrame:
    """Unified entry point; always emits column named SIGNAL for downstream modules."""
    if bin_mode:
        sig = compute_mean_bin(reads).rename({MEAN_BIN: SIGNAL})
    else:
        sig = compute_frequency_signal(reads)
    return sig
