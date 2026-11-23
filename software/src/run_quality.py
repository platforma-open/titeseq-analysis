"""Counts of fit outcomes and of failure reasons, for the block's run-quality graphs.

Counted on the lowest grain the run fitted: one row per variant, and per replicate when the run has
replicates. A combined value's failure reason only says that its parts failed.

The graphs count the variants the experiment set out to measure. With a known-variant list from
the upstream block, those are the variants that match a known one exactly. Without one, they are
the variants with the reads and concentration points for a fit to be attempted: on raw reads the
rest are mostly sequencing errors and outnumber every class. The raw counts keep every variant.
"""

from __future__ import annotations

import polars as pl

from constants import COL_CLONOTYPE, COL_REPLICATE

COL_REASON = "failureReason"
COL_COUNT = "count"
# The failure reasons that mean a fit was never attempted.
NOT_ATTEMPTED = ("insufficient_reads", "insufficient_points")
# The replicate value of the counts over all replicates together, unless a replicate has that name.
ALL_REPLICATES = "All replicates"
SCOPE_KNOWN = "known"
SCOPE_ATTEMPTED = "attempted"


def fit_grain(outputs: dict) -> pl.DataFrame:
    """The frame holding one row per fit the run made."""
    for name in ("per_replicate_nt", "per_replicate"):
        if name in outputs:
            return outputs[name]
    return outputs["per_clonotype"]


def combined_label(frame: pl.DataFrame) -> str:
    """The replicate value for the counts over all replicates: one no replicate of `frame` has."""
    label = ALL_REPLICATES
    if COL_REPLICATE not in frame.columns:
        return label
    taken = set(frame[COL_REPLICATE].cast(pl.Utf8).to_list())
    while label in taken:
        label += " (combined)"
    return label


def _count_by(frame: pl.DataFrame, columns: list[str], combined: str) -> pl.DataFrame:
    """Rows of `frame` per value of `columns`, per replicate and under `combined` for all replicates."""
    count = pl.len().cast(pl.Int64).alias(COL_COUNT)
    if COL_REPLICATE not in frame.columns:
        return frame.group_by(columns).agg(count).sort(columns, nulls_last=True)
    keys = [COL_REPLICATE, *columns]
    each = frame.with_columns(pl.col(COL_REPLICATE).cast(pl.Utf8)).group_by(keys).agg(count)
    together = frame.group_by(columns).agg(count).with_columns(pl.lit(combined).alias(COL_REPLICATE))
    together = together.sort(columns, nulls_last=True).select([*keys, COL_COUNT])
    return pl.concat([together, each.sort(keys, nulls_last=True)])


def _counted(frame: pl.DataFrame, column: str, combined: str) -> pl.DataFrame:
    return _count_by(frame, [column], combined)


def read_known(path: str | None) -> set[str] | None:
    """The variants that match a known one exactly, from a workflow-built TSV [key, distance].

    None when the run sent no table: the upstream block carries no known-variant list.
    """
    if path is None:
        return None
    frame = pl.read_csv(path, separator="\t", infer_schema=False)
    exact = frame.filter(pl.col("distance").cast(pl.Int64, strict=False) == 0)
    return set(exact["key"].to_list())


def measured(frame: pl.DataFrame, known: set[str] | None) -> pl.DataFrame:
    """The fits the graphs count."""
    if known is not None:
        return frame.filter(pl.col(COL_CLONOTYPE).cast(pl.Utf8).is_in(list(known)))
    return frame.filter(~pl.col("fitFailureReason").fill_null("").is_in(NOT_ATTEMPTED))


def _labelled(frame: pl.DataFrame) -> pl.DataFrame:
    return frame.with_columns(
        pl.col("fitFailureReason").map_elements(reason_label, return_dtype=pl.Utf8).alias(COL_REASON)
    )


def raw_counts(frame: pl.DataFrame, combined: str = ALL_REPLICATES) -> pl.DataFrame:
    """[replicate?, affinityClass, failureReason, count] over every fit; the reason is null unless it failed."""
    return _count_by(_labelled(frame), ["affinityClass", COL_REASON], combined)


def outcome_counts(frame: pl.DataFrame, combined: str = ALL_REPLICATES) -> pl.DataFrame:
    """[replicate?, affinityClass, count]."""
    return _counted(frame, "affinityClass", combined)


def reason_label(code: str) -> str:
    """A failure reason as a legend reads it: "insufficient_points" -> "Insufficient points"."""
    return code.replace("_", " ").capitalize()


def failure_reason_counts(frame: pl.DataFrame, combined: str = ALL_REPLICATES) -> pl.DataFrame:
    """[replicate?, failureReason, count] over the failed fits."""
    failed = frame.filter(pl.col("affinityClass") == "Failed").with_columns(
        pl.col("fitFailureReason")
        .fill_null("unknown")
        .map_elements(reason_label, return_dtype=pl.Utf8)
        .alias(COL_REASON)
    )
    return _counted(failed, COL_REASON, combined)


def summary(
    outcomes: pl.DataFrame,
    reasons: pl.DataFrame,
    raw: pl.DataFrame,
    known: set[str] | None,
    combined: str = ALL_REPLICATES,
) -> dict:
    """The tables as JSON-ready rows, and what the graphs count."""
    return {
        "scope": SCOPE_ATTEMPTED if known is None else SCOPE_KNOWN,
        "allReplicates": combined,
        "leftOutReasons": [reason_label(r) for r in NOT_ATTEMPTED] if known is None else [],
        "outcomes": outcomes.to_dicts(),
        "failureReasons": reasons.to_dicts(),
        "raw": raw.to_dicts(),
    }
