"""Sample replicates combined into one value and one class per variant.

Each replicate is a titration of its own: its samples are fitted on their own reads, its gate depths
and sort fractions are its own, and it produces its own K_D per variant. This module takes those
per-replicate results and reduces them to one.

  entered      a replicate with a usable K_D: measured, or censored at a wall of the titration
  measured     the fit put K_D inside the titrated range
  below/above  K_D ran past a wall; the replicate enters at that wall, as Adams's and Phillips's
               code both do
  failed       the replicate produced no K_D
  filtered     the replicate measured a K_D that a pre-combine filter rejected; left out of the
               mean and counted. Every filter is off by default
  value        the plain mean of the entered replicates' log10 K_D
  SE           by `error_model`:
                 sd_over_sqrt_n  sample SD (ddof 1) / sqrt(n); for n = 1 the replicate's own fit
                                 SE, where the equation reports one, else empty
                 phillips        population SD (ddof 0) / sqrt(n + 1), which is 0 for n = 1.
                                 Reproduces the SEM published with Phillips 2021
  SD           sample SD (ddof 1) of the entered replicates, whatever the error model; empty for
               n = 1
  bound        no entered replicate measured a K_D. One censoring side must hold at least
               `censored_majority` of them and more than the other side; the value is that side's
               wall. Otherwise the sides conflict and the variant fails

The class is decided once, on the combined value, exactly as the synonymous combination decides it:
Failed with no entered replicate, Partial for a bound, Good for an SE at or below KD_SD_GOOD,
Partial above it, and the single replicate's own class when n = 1 and there is no SE.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from typing import Literal

import polars as pl

from classify import Classification
from constants import (
    COL_CLONOTYPE,
    COL_REPLICATE,
    DEFAULT_CENSORED_MAJORITY,
    KD_SD_GOOD,
    ReplicateErrorModel,
)
from output_build import PER_CLONOTYPE_SCHEMA
from synonymous import COL_AA, PER_VARIANT_SCHEMA, copy_from_fit

ReplicateState = Literal["measured", "below", "above", "failed", "filtered"]
STATES: tuple[ReplicateState, ...] = ("measured", "below", "above", "failed", "filtered")
ENTERED: tuple[ReplicateState, ...] = ("measured", "below", "above")

REASON_ALL_FILTERED = "all_replicates_filtered"
REASON_CONFLICTING = "conflicting_censoring"


@dataclass(frozen=True)
class Filters:
    """The pre-combine filters, each off when None. Both read a fit, so both apply to a measured
    replicate only: a censored value is set by the wall, not by the fit."""

    fit_error_cutoff: float | None = None
    r2_floor: float | None = None

    @property
    def active(self) -> bool:
        return self.fit_error_cutoff is not None or self.r2_floor is not None

    def rejects(self, log10_kd_se: float | None, r2: float | None) -> bool:
        if self.fit_error_cutoff is not None:
            if log10_kd_se is None or not math.isfinite(log10_kd_se) or log10_kd_se > self.fit_error_cutoff:
                return True
        if self.r2_floor is not None:
            if r2 is None or not math.isfinite(r2) or r2 < self.r2_floor:
                return True
        return False


NO_FILTERS = Filters()


@dataclass(frozen=True)
class ReplicateValue:
    """One replicate's reported K_D for one variant, as it enters the mean."""

    replicate: str
    state: ReplicateState
    log10_kd: float | None
    log10_kd_se: float | None
    affinity_class: str
    r2: float | None = None


@dataclass(frozen=True)
class Merged:
    log10_kd: float | None
    log10_kd_se: float | None
    log10_kd_sd: float | None
    classification: Classification
    out_of_range: bool | None
    counts: dict[ReplicateState, int]


def replicate_value(
    replicate: str,
    row: dict,
    *,
    log10_lo: float | None,
    log10_hi: float | None,
    filters: Filters = NO_FILTERS,
) -> ReplicateValue:
    """One replicate's fit row as it enters the mean. `row` is a per-clonotype or per-variant row."""
    copy = copy_from_fit(row, log10_lo, log10_hi)
    r2 = row.get("r2")
    state: ReplicateState = copy.state
    if state == "measured" and filters.rejects(copy.log10_kd_se, r2):
        state = "filtered"
    return ReplicateValue(replicate, state, copy.log10_kd, copy.log10_kd_se, copy.affinity_class, r2)


def _spread_se(values: list[float], error_model: ReplicateErrorModel) -> float:
    n = len(values)
    if error_model == "phillips":
        return statistics.pstdev(values) / math.sqrt(n + 1)
    return statistics.stdev(values) / math.sqrt(n)


def _bound(entered: list[ReplicateValue], counts: dict[ReplicateState, int], majority: float) -> Merged:
    """No entered replicate measured a K_D: the value is a wall, or the sides conflict.

    Averaging opposite walls would land inside the titrated range and mean nothing, so one side
    must hold the majority. The threshold is the block's own rule; no paper defines one.
    """
    n = len(entered)
    share = {side: sum(1 for v in entered if v.state == side) / n for side in ("below", "above")}
    for side, other in (("below", "above"), ("above", "below")):
        if share[side] >= majority and share[side] > share[other]:
            wall = next(v.log10_kd for v in entered if v.state == side)
            return Merged(wall, None, None, Classification("Partial", None), True, counts)
    return Merged(None, None, None, Classification("Failed", REASON_CONFLICTING), None, counts)


def combine_replicates(
    values: list[ReplicateValue],
    *,
    error_model: ReplicateErrorModel = "sd_over_sqrt_n",
    censored_majority: float = DEFAULT_CENSORED_MAJORITY,
) -> Merged:
    """One variant's value over its replicates, and its class."""
    counts: dict[ReplicateState, int] = dict.fromkeys(STATES, 0)
    for v in values:
        counts[v.state] += 1
    entered = [v for v in values if v.state in ENTERED and v.log10_kd is not None]
    n = len(entered)
    if n == 0:
        reason = REASON_ALL_FILTERED if counts["filtered"] > 0 else None
        return Merged(None, None, None, Classification("Failed", reason), None, counts)

    if not any(v.state == "measured" for v in entered):
        return _bound(entered, counts, censored_majority)

    points = [v.log10_kd for v in entered]
    mean = math.fsum(points) / n
    if n == 1:
        se = 0.0 if error_model == "phillips" else entered[0].log10_kd_se
        sd = None
    else:
        se = _spread_se(points, error_model)
        sd = statistics.stdev(points)
    if se is not None:
        cls = "Good" if math.isfinite(se) and se <= KD_SD_GOOD else "Partial"
        return Merged(mean, se, sd, Classification(cls, None), False, counts)
    return Merged(mean, None, sd, Classification(entered[0].affinity_class, None), False, counts)


@dataclass(frozen=True)
class ReplicateInputs:
    """The replicate column and every setting the combination reads."""

    column: str = COL_REPLICATE
    error_model: ReplicateErrorModel = "sd_over_sqrt_n"
    filters: Filters = field(default_factory=Filters)
    censored_majority: float = DEFAULT_CENSORED_MAJORITY


# How each replicate's state reads in the per-replicate output.
STATE_LABELS: dict[ReplicateState, str] = {
    "measured": "measured",
    "below": "censoredBelow",
    "above": "censoredAbove",
    "failed": "failed",
    "filtered": "filtered",
}

PER_REPLICATE_COLUMNS: dict[str, pl.DataType] = {
    COL_REPLICATE: pl.Utf8,
    "kd": pl.Float64,
    "kdLog10Se": pl.Float64,
    "r2": pl.Float64,
    "affinityClass": pl.Utf8,
    "kdState": pl.Utf8,
}

# Counted per variant: how many replicates reached each state. `replicateCount` is the number that
# entered the mean, which is the measured and censored ones together.
COUNT_SCHEMA: dict[str, pl.DataType] = {
    "replicateCount": pl.Int64,
    "replicateMeasuredCount": pl.Int64,
    "replicateCensoredBelowCount": pl.Int64,
    "replicateCensoredAboveCount": pl.Int64,
    "replicateFailedCount": pl.Int64,
    "replicateFilteredCount": pl.Int64,
}

Collected = dict[str, list[tuple[ReplicateValue, dict]]]


def collect_values(
    frames: list[tuple[str, pl.DataFrame]],
    key_col: str,
    *,
    log10_lo: float | None,
    log10_hi: float | None,
    filters: Filters = NO_FILTERS,
) -> Collected:
    """Each variant's per-replicate values, paired with the fit row each came from."""
    collected: Collected = {}
    for label, frame in frames:
        for row in frame.with_columns(pl.col(key_col).cast(pl.Utf8)).iter_rows(named=True):
            value = replicate_value(label, row, log10_lo=log10_lo, log10_hi=log10_hi, filters=filters)
            collected.setdefault(row[key_col], []).append((value, row))
    return collected


def merge_all(
    collected: Collected,
    *,
    error_model: ReplicateErrorModel = "sd_over_sqrt_n",
    censored_majority: float = DEFAULT_CENSORED_MAJORITY,
) -> dict[str, Merged]:
    """One `Merged` per variant, keys in sorted order."""
    return {
        key: combine_replicates(
            [v for v, _ in collected[key]],
            error_model=error_model,
            censored_majority=censored_majority,
        )
        for key in sorted(collected)
    }


def build_per_replicate_frame(collected: Collected, key_col: str) -> pl.DataFrame:
    """One row per (variant, replicate): the value that replicate contributed, and its state.

    `kd` is the value as it enters the mean, so a censored replicate reports the wall it passed
    and a filtered one reports the K_D the filter rejected.
    """
    schema: dict[str, pl.DataType] = {key_col: pl.Utf8, **PER_REPLICATE_COLUMNS}
    cols: dict[str, list] = {name: [] for name in schema}
    for key in sorted(collected):
        for value, _row in collected[key]:
            cols[key_col].append(key)
            cols[COL_REPLICATE].append(value.replicate)
            cols["kd"].append(None if value.log10_kd is None else 10.0**value.log10_kd)
            cols["kdLog10Se"].append(value.log10_kd_se)
            cols["r2"].append(value.r2)
            cols["affinityClass"].append(value.affinity_class)
            cols["kdState"].append(STATE_LABELS[value.state])
    return pl.DataFrame(cols, schema=schema)


def _counts_of(merged: Merged) -> list[int]:
    c = merged.counts
    entered = c["measured"] + c["below"] + c["above"]
    return [entered, c["measured"], c["below"], c["above"], c["failed"], c["filtered"]]


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


def build_combined_clonotype_frame(collected: Collected, merged: dict[str, Merged]) -> pl.DataFrame:
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
        cols["kdLog10Sd"].append(m.log10_kd_sd)
        for name, count in zip(COUNT_SCHEMA, _counts_of(m)):
            cols[name].append(count)
    return pl.DataFrame(cols, schema=COMBINED_CLONOTYPE_SCHEMA)


# The counts the synonymous combination produced inside each replicate, summed over them.
_NT_COUNT_FIELDS = ("ntVariants", "ntMeasured", "ntCensoredBelow", "ntCensoredAbove", "ntFailed")

COMBINED_VARIANT_SCHEMA: dict[str, pl.DataType] = {**PER_VARIANT_SCHEMA, **COUNT_SCHEMA}


def build_combined_variant_frame(collected: Collected, merged: dict[str, Merged]) -> pl.DataFrame:
    """The amino-acid frame, combined over the replicates after the synonymous combination.

    The nucleotide counts are the sum over the replicates: each replicate fitted its own copies.
    """
    cols: dict[str, list] = {name: [] for name in COMBINED_VARIANT_SCHEMA}
    for key, m in merged.items():
        rows = collected[key]
        cols[COL_AA].append(key)
        cols["kd"].append(None if m.log10_kd is None else 10.0**m.log10_kd)
        cols["kdLog10Se"].append(m.log10_kd_se)
        cols["kdLog10Sd"].append(m.log10_kd_sd)
        cols["affinityClass"].append(m.classification.affinity_class)
        cols["kdOutOfRange"].append(m.out_of_range)
        for name in _NT_COUNT_FIELDS:
            cols[name].append(sum(int(row[name] or 0) for _, row in rows))
        for name, count in zip(COUNT_SCHEMA, _counts_of(m)):
            cols[name].append(count)
    return pl.DataFrame(cols, schema=COMBINED_VARIANT_SCHEMA)
