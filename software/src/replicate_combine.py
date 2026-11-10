"""Sample replicates combined into one value and one class per variant.

Each replicate is a titration of its own: its samples are fitted on their own reads, its gate depths
and sort fractions are its own, and it produces its own K_D per variant. This module takes those
per-replicate results and reduces them to one.

  entered      a replicate with a usable K_D: measured, or censored at a wall of the titration
  measured     the fit put K_D inside the titrated range
  below/above  K_D ran past a wall; the replicate enters at that wall, as Adams's and Phillips's
               code both do
  failed       the replicate produced no K_D
  value        the plain mean of the entered replicates' log10 K_D
  SE           for n >= 2, max(spread, fit-error term), as for synonymous variants
               (`classify.rate_mean`). The spread is the sample SD (ddof 1) / sqrt(n). The
               fit-error term is the RMS of the SEs the replicates carry, over sqrt(n), n counting
               every entered replicate; with no SE the SE is the spread. For n = 1, the replicate's
               own SE where it carries one, else empty
  SD           sample SD (ddof 1) of the entered replicates; empty for n = 1
  bound        no entered replicate measured a K_D. One censoring side must hold at least half of
               them and more than the other side; the value is that side's wall. Otherwise the
               sides conflict and the variant fails

The class is decided once, on the combined value: Failed with no entered replicate; for a bound,
the wall it sits at under "log_fluorescence", and under "hill_gate_index" the class its replicates
share, else Partial (`classify.bound_class`); Good for an SE at or below KD_SD_GOOD with more than
half the entered replicates measured and, under "hill_gate_index" or when no replicate carries an
SE, more than half of them rated Good; Partial otherwise; and the single replicate's own class when
n = 1 and there is no SE.

The R² gate (`FitParams.r2_gate_min`) runs inside each replicate's own fit and rates the fit
Failed, so the replicate enters as "failed".
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, replace
from functools import partial

import polars as pl

from classify import Classification, Combined, Member, MemberState, combine_members, member_from_fit
from constants import COL_CLONOTYPE, COL_REPLICATE
from output_build import PER_CLONOTYPE_SCHEMA
from synonymous import COL_AA, PER_VARIANT_SCHEMA

ENTERED: tuple[MemberState, ...] = ("measured", "below", "above")

REASON_ALL_FAILED = "all_replicates_failed"

# The spread the SE reads, as the run summary names it.
ERROR_MODEL = "sd_over_sqrt_n"


@dataclass(frozen=True, kw_only=True)
class ReplicateValue(Member):
    """One replicate's reported K_D for one variant, as it enters the mean."""

    replicate: str
    r2: float | None = None
    # The replicate's own fit failure reason; None unless its fit rated Failed.
    failure_reason: str | None = None


def replicate_value(replicate: str, row: dict, *, log10_lo: float | None, log10_hi: float | None) -> ReplicateValue:
    """One replicate's fit row as it enters the mean. `row` is a per-clonotype row."""
    fit = member_from_fit(row, log10_lo, log10_hi)
    return ReplicateValue(
        fit.state,
        fit.log10_kd,
        fit.log10_kd_se,
        fit.affinity_class,
        replicate=replicate,
        r2=row.get("r2"),
        failure_reason=row.get("fitFailureReason"),
    )


# One variant's value over its replicates, and its class (`classify.combine_members`).
combine_replicates = partial(combine_members, all_failed_reason=REASON_ALL_FAILED)


@dataclass(frozen=True)
class ReplicateInputs:
    """The column that names each sample's replicate."""

    column: str = COL_REPLICATE


PER_REPLICATE_COLUMNS: dict[str, pl.DataType] = {
    COL_REPLICATE: pl.Utf8,
    "replicateReads": pl.Int64,
    "kd": pl.Float64,
    "kdNegLog10": pl.Float64,
    "kdLog10Se": pl.Float64,
    "r2": pl.Float64,
    "affinityClass": pl.Utf8,
    "fitFailureReason": pl.Utf8,
}

# Counted per variant: how many replicates reached each state. `replicateCount` is the number that
# entered the mean, which is the measured and censored ones together.
COUNT_SCHEMA: dict[str, pl.DataType] = {
    "replicateCount": pl.Int64,
    "replicateMeasuredCount": pl.Int64,
    "replicateCensoredBelowCount": pl.Int64,
    "replicateCensoredAboveCount": pl.Int64,
    "replicateFailedCount": pl.Int64,
}

Collected = dict[str, list[tuple[ReplicateValue, dict]]]


def collect_values(
    frames: list[tuple[str, pl.DataFrame]],
    key_col: str,
    *,
    log10_lo: float | None,
    log10_hi: float | None,
) -> Collected:
    """Each variant's per-replicate values, paired with the fit row each came from."""
    collected: Collected = {}
    for label, frame in frames:
        for row in frame.with_columns(pl.col(key_col).cast(pl.Utf8)).iter_rows(named=True):
            value = replicate_value(label, row, log10_lo=log10_lo, log10_hi=log10_hi)
            collected.setdefault(row[key_col], []).append((value, row))
    return collected


def merge_all(collected: Collected, *, wall_classes: bool = True, fits_report_se: bool = True) -> dict[str, Combined]:
    """One combined value per variant, keys in sorted order."""
    return {
        key: combine_replicates(
            [v for v, _ in collected[key]], wall_classes=wall_classes, fits_report_se=fits_report_se
        )
        for key in sorted(collected)
    }


def build_per_replicate_frame(
    collected: Collected,
    key_col: str,
    reads: dict[tuple[str, str], int] | None = None,
    extra: Mapping[str, pl.DataType] | None = None,
) -> pl.DataFrame:
    """One row per (variant, replicate): the value that replicate contributed, and its state.

    `kd` is the value as it enters the mean, so a censored replicate reports the wall it passed.
    `reads` maps (variant, replicate) to every read of the variant in that replicate, which is what
    separates a wide value from a wrong one.
    `extra` names further columns copied as they are from each value's source row.
    """
    extra = extra or {}
    schema: dict[str, pl.DataType] = {key_col: pl.Utf8, **PER_REPLICATE_COLUMNS, **extra}
    cols: dict[str, list] = {name: [] for name in schema}
    for key in sorted(collected):
        for value, row in collected[key]:
            for name in extra:
                cols[name].append(row[name])
            cols[key_col].append(key)
            cols[COL_REPLICATE].append(value.replicate)
            cols["replicateReads"].append(None if reads is None else reads.get((key, value.replicate)))
            cols["kd"].append(None if value.log10_kd is None else 10.0**value.log10_kd)
            # The published form, as the combined value reports it: larger is tighter.
            cols["kdNegLog10"].append(None if value.log10_kd is None else -value.log10_kd)
            cols["kdLog10Se"].append(value.log10_kd_se)
            cols["r2"].append(value.r2)
            cols["affinityClass"].append(value.affinity_class)
            cols["fitFailureReason"].append(value.failure_reason)
    return pl.DataFrame(cols, schema=schema)


def _counts_of(merged: Combined) -> list[int]:
    c = merged.counts
    entered = c["measured"] + c["below"] + c["above"]
    return [entered, c["measured"], c["below"], c["above"], c["failed"]]


def _mean_over_entered(rows: list[tuple[ReplicateValue, dict]], field_name: str) -> float | None:
    """The mean of `field_name` over the replicates that entered the mean; None when none has it."""
    values = [
        row[field_name]
        for value, row in rows
        if value.state in ENTERED and row.get(field_name) is not None and math.isfinite(row[field_name])
    ]
    return math.fsum(values) / len(values) if values else None


COMBINED_CLONOTYPE_SCHEMA: dict[str, pl.DataType] = {
    **PER_CLONOTYPE_SCHEMA,
    "kdLog10Sd": pl.Float64,
    **COUNT_SCHEMA,
}


def build_combined_clonotype_frame(collected: Collected, merged: dict[str, Combined]) -> pl.DataFrame:
    """The per-clonotype frame the rest of the pipeline reads, combined over the replicates.

    `hillCoefficient` and `r2` are the mean over the replicates that entered: no equation defines
    a combined fit, and the mean keeps both columns readable as a fit-quality summary.
    """
    cols: dict[str, list] = {name: [] for name in COMBINED_CLONOTYPE_SCHEMA}
    for key, m in merged.items():
        rows = collected[key]
        cols[COL_CLONOTYPE].append(key)
        cols["kd"].append(None if m.log10_kd is None else 10.0**m.log10_kd)
        cols["hillCoefficient"].append(_mean_over_entered(rows, "hillCoefficient"))
        cols["r2"].append(_mean_over_entered(rows, "r2"))
        cols["affinityClass"].append(m.classification.affinity_class)
        cols["fitFailureReason"].append(m.classification.failure_reason)
        cols["kdOutOfRange"].append(m.out_of_range)
        cols["kdLog10Se"].append(m.log10_kd_se)
        cols["kdRangeSide"].append(m.range_side)
        cols["kdLog10Sd"].append(m.log10_kd_sd)
        for name, count in zip(COUNT_SCHEMA, _counts_of(m)):
            cols[name].append(count)
    return pl.DataFrame(cols, schema=COMBINED_CLONOTYPE_SCHEMA)


# The counts the synonymous combination produced inside each replicate, summed over them.
_NT_STATE_COUNT_FIELDS = ("ntMeasured", "ntCensoredBelow", "ntCensoredAbove", "ntFailed")

COMBINED_VARIANT_SCHEMA: dict[str, pl.DataType] = {**PER_VARIANT_SCHEMA, **COUNT_SCHEMA}


def _has_measured_majority(row: dict) -> bool:
    """True when more than half of one replicate's synonymous variants in its mean were measured."""
    measured = int(row["ntMeasured"] or 0)
    at_a_wall = int(row["ntCensoredBelow"] or 0) + int(row["ntCensoredAbove"] or 0)
    return measured > at_a_wall


def require_synonymous_majority(collected: Collected, merged: dict[str, Combined]) -> dict[str, Combined]:
    """Good amino-acid values downgraded to Partial unless, in more than half of the replicates that
    entered the mean, more than half of the synonymous variants' fits were measured.

    A replicate's amino-acid value with a measured synonymous variant enters the replicates' mean as
    measured, however many of its synonymous variants sat at a wall. Two such replicates agree on the
    walls and shrink the SE, so the majority rule is applied again on the synonymous variants. It is
    counted inside each replicate: pooled over the replicates, one with many measured synonymous
    variants would carry the others. A bound is left as it is: it has no SE for the rule to guard,
    and it keeps the class it has without replicates.
    """
    out: dict[str, Combined] = {}
    for key, m in merged.items():
        if m.classification.affinity_class == "Good" and not m.out_of_range:
            rows = [row for value, row in collected[key] if value.state in ENTERED]
            if 2 * sum(_has_measured_majority(row) for row in rows) <= len(rows):
                m = replace(m, classification=Classification("Partial", None))
        out[key] = m
    return out


def build_combined_variant_frame(
    collected: Collected, merged: dict[str, Combined], nt_variants: dict[str, int]
) -> pl.DataFrame:
    """The amino-acid frame, combined over the replicates after the synonymous combination.

    `nt_variants` is each amino-acid variant's count of distinct synonymous variants fitted in any
    replicate. The counts by state are the sum over the replicates: each replicate fitted its own
    synonymous variants.
    """
    cols: dict[str, list] = {name: [] for name in COMBINED_VARIANT_SCHEMA}
    for key, m in merged.items():
        rows = collected[key]
        cols[COL_AA].append(key)
        cols["kd"].append(None if m.log10_kd is None else 10.0**m.log10_kd)
        cols["kdNegLog10"].append(None if m.log10_kd is None else -m.log10_kd)
        cols["kdLog10Se"].append(m.log10_kd_se)
        cols["kdLog10Sd"].append(m.log10_kd_sd)
        cols["affinityClass"].append(m.classification.affinity_class)
        cols["fitFailureReason"].append(m.classification.failure_reason)
        cols["kdOutOfRange"].append(m.out_of_range)
        cols["ntVariants"].append(nt_variants.get(key, 0))
        for name in _NT_STATE_COUNT_FIELDS:
            cols[name].append(sum(int(row[name] or 0) for _, row in rows))
        for name, count in zip(COUNT_SCHEMA, _counts_of(m)):
            cols[name].append(count)
    return pl.DataFrame(cols, schema=COMBINED_VARIANT_SCHEMA)
