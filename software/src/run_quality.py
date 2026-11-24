"""Counts of fit outcomes and of failure reasons, for the block's run-quality graphs.

Counted on the lowest grain the run fitted: one row per variant, and per replicate when the run has
replicates. A combined value's failure reason only says that its parts failed.

The graphs count the variants the experiment set out to measure. With a known-variant list from
the upstream block, those are the variants that match a known one exactly. Without one, they are
the variants with the reads and concentration points for a fit to be attempted: on raw reads the
rest are mostly sequencing errors and outnumber every class. The raw counts keep every variant.

The two TSV tables hold each replicate's counts. The JSON summary holds those and the counts over
all replicates together.
"""

from __future__ import annotations

import json

import polars as pl

from constants import COL_CLONOTYPE, COL_REPLICATE

COL_CODE = "fitFailureReason"
COL_REASON = "failureReason"
COL_COUNT = "count"
# The failure reasons of a variant without the reads or the concentration points for a fit.
NOT_ATTEMPTED = ("insufficient_reads", "insufficient_points")
# The labels that do not follow from the code. The workflow and the Run Quality page list every
# label a fit can report: keep the three in step.
REASON_LABELS = {
    "low_r2": "Low R²",
    "high_kd_se": "High Kd,app SE",
    "n_out_of_range": "Hill coefficient out of range",
}
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
    return frame.filter(~pl.col(COL_CODE).fill_null("").is_in(NOT_ATTEMPTED))


def reason_label(code: str) -> str:
    """A failure reason as a legend reads it: "insufficient_points" -> "Insufficient points"."""
    return REASON_LABELS.get(code) or code.replace("_", " ").capitalize()


def _count_by_reason(frame: pl.DataFrame, columns: list[str], combined: str) -> pl.DataFrame:
    """`_count_by` over `columns` and the failure reason, the counted reasons then labelled."""
    counts = _count_by(frame, [*columns, COL_CODE], combined).rename({COL_CODE: COL_REASON})
    return counts.with_columns(pl.col(COL_REASON).cast(pl.Utf8).map_elements(reason_label, return_dtype=pl.Utf8))


def raw_counts(frame: pl.DataFrame, combined: str = ALL_REPLICATES) -> pl.DataFrame:
    """[replicate?, affinityClass, failureReason, count] over every fit; the reason is null unless it failed."""
    return _count_by_reason(frame, ["affinityClass"], combined)


def outcome_counts(frame: pl.DataFrame, combined: str = ALL_REPLICATES) -> pl.DataFrame:
    """[replicate?, affinityClass, count]."""
    return _count_by(frame, ["affinityClass"], combined)


def failure_reason_counts(frame: pl.DataFrame, combined: str = ALL_REPLICATES) -> pl.DataFrame:
    """[replicate?, failureReason, count] over the failed fits. Every failed fit carries a reason."""
    return _count_by_reason(frame.filter(pl.col("affinityClass") == "Failed"), [], combined)


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


def per_replicate(counts: pl.DataFrame, combined: str) -> pl.DataFrame:
    """`counts` without the rows over all replicates."""
    if COL_REPLICATE not in counts.columns:
        return counts
    return counts.filter(pl.col(COL_REPLICATE) != combined)


def write_run_quality(outputs: dict, known_path: str | None, prefix: str) -> None:
    """Write <prefix>outcomes.tsv, <prefix>failure_reasons.tsv and <prefix>summary.json.

    The tables hold each replicate's counts; the summary also holds the counts over all replicates.
    """
    grain = fit_grain(outputs)
    known = read_known(known_path)
    counted = measured(grain, known)
    combined = combined_label(grain)
    outcomes, reasons = outcome_counts(counted, combined), failure_reason_counts(counted, combined)
    per_replicate(outcomes, combined).write_csv(prefix + "outcomes.tsv", separator="\t")
    per_replicate(reasons, combined).write_csv(prefix + "failure_reasons.tsv", separator="\t")
    with open(prefix + "summary.json", "w") as f:
        json.dump(summary(outcomes, reasons, raw_counts(grain, combined), known, combined), f)
