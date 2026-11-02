"""Synonymous nucleotide copies combined into one value per amino-acid variant, and its parent delta.

Every nucleotide variant is fitted on its own before this module runs. Only the per-copy results are
combined here; raw counts of two copies are never pooled.

The combination is the arithmetic mean of the copies' log10 K_D (Starr 2020, 2022):

  copies       the fitted nucleotide variants that map to exactly one amino-acid variant
  failed       a copy rated Failed, or with no usable K_D; left out of the mean and counted
  censored     a copy whose K_D sits outside the titrated range; enters at the wall it passed, and is
               counted by direction ("below" the lowest concentration, "above" the highest)
  measured     every other copy; enters at its own log10 K_D
  SE           for n >= 2 copies in the mean, the larger of two terms: the spread, sample SD / sqrt(n),
               and the fit-error term, sqrt(sum of se_i^2) / n over the copies that have a fit SE
               (null when none has one, and the SE is then the spread); for n = 1, the copy's own
               fit SE where the equation reports one, else null
  bound        every copy in the mean censored on the same side; the value is that wall, with no SE

The class is decided once, on the combined value:

  Failed       no copy in the mean
  Partial      a bound
  Good         SE at or below KD_SD_GOOD
  Partial      SE above KD_SD_GOOD, or not finite
  copy's own   n = 1 and the equation reports no SE: the combined value is that copy's fit, rated by
               its equation's own rule

The parent delta is log10 K_D of the variant minus that of its parent, both at amino-acid grain, with
error sqrt(se_v^2 + se_p^2). A parent's own row is 0 with error 0. A parent is the one amino-acid variant
with zero mutations linked to it; a parent with none or several cannot be identified.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass
from typing import Literal

import polars as pl

from classify import Classification
from constants import COL_CLONOTYPE, KD_SD_GOOD

COL_NT = "ntKey"
COL_AA = "aaKey"
COL_PARENT = "parentId"
COL_MUTATIONS = "mutationCount"

REASON_NO_UNMUTATED = "no-unmutated-variant"
REASON_SEVERAL_UNMUTATED = "several-unmutated-variants"

CopyState = Literal["measured", "below", "above", "failed"]

PER_VARIANT_SCHEMA: dict[str, pl.DataType] = {
    COL_AA: pl.Utf8,
    "kd": pl.Float64,
    "kdLog10Se": pl.Float64,
    "affinityClass": pl.Utf8,
    "kdOutOfRange": pl.Boolean,
    "ntVariants": pl.Int64,
    "ntMeasured": pl.Int64,
    "ntCensoredBelow": pl.Int64,
    "ntCensoredAbove": pl.Int64,
    "ntFailed": pl.Int64,
}
DELTA_SCHEMA: dict[str, pl.DataType] = {
    "kdLog10ParentDelta": pl.Float64,
    "kdLog10ParentDeltaSe": pl.Float64,
}


@dataclass(frozen=True)
class SynonymousInputs:
    """The nt -> aa linker, and the amino-acid parent columns when both resolve."""

    nt_to_aa: pl.DataFrame  # [ntKey, aaKey]
    aa_mutations: pl.DataFrame | None = None  # [aaKey, mutationCount]
    aa_parent_links: pl.DataFrame | None = None  # [aaKey, parentId]

    def __post_init__(self) -> None:
        if (self.aa_mutations is None) != (self.aa_parent_links is None):
            raise ValueError("the amino-acid mutation counts and parent links come together or not at all")


@dataclass(frozen=True)
class Copy:
    state: CopyState
    log10_kd: float | None
    log10_kd_se: float | None
    affinity_class: str


@dataclass(frozen=True)
class Combined:
    log10_kd: float | None
    log10_kd_se: float | None
    classification: Classification
    out_of_range: bool | None
    counts: dict[CopyState, int]


def _keys(frame: pl.DataFrame, *cols: str) -> pl.DataFrame:
    return frame.select([pl.col(c).cast(pl.Utf8) for c in cols]).drop_nulls().unique()


def resolve_nt_to_aa(nt_to_aa: pl.DataFrame, fitted: list[str]) -> tuple[dict[str, str], int, int]:
    """Each fitted nt variant's single aa variant; then how many map to none and to several."""
    per_nt = _keys(nt_to_aa, COL_NT, COL_AA).group_by(COL_NT).agg(pl.col(COL_AA).unique())
    targets = dict(zip(per_nt[COL_NT].to_list(), per_nt[COL_AA].to_list()))
    mapping: dict[str, str] = {}
    unmapped = ambiguous = 0
    for nt in fitted:
        aas = targets.get(nt)
        if aas is None:
            unmapped += 1
        elif len(aas) > 1:
            ambiguous += 1
        else:
            mapping[nt] = aas[0]
    return mapping, unmapped, ambiguous


def copy_from_fit(row: dict, log10_lo: float | None, log10_hi: float | None) -> Copy:
    """One nt variant's fit as it enters the mean. `row` is a per-clonotype output row."""
    kd = row["kd"]
    cls = row["affinityClass"]
    if cls == "Failed" or kd is None or not math.isfinite(kd) or kd <= 0:
        return Copy("failed", None, None, cls)
    value = math.log10(kd)
    if row["kdOutOfRange"] and log10_lo is not None and log10_hi is not None:
        if value <= (log10_lo + log10_hi) / 2:
            return Copy("below", log10_lo, None, cls)
        return Copy("above", log10_hi, None, cls)
    se = row["kdLog10Se"]
    return Copy("measured", value, se if se is not None and math.isfinite(se) else None, cls)


def _floored_se(entered: list[Copy]) -> float:
    """max(spread, fit-error term) over n >= 2 copies; the spread alone when no copy has a fit SE."""
    n = len(entered)
    spread = statistics.stdev([c.log10_kd for c in entered]) / math.sqrt(n)
    fit_ses = [c.log10_kd_se for c in entered if c.log10_kd_se is not None]
    if not fit_ses:
        return spread
    return max(spread, math.sqrt(math.fsum(se * se for se in fit_ses)) / n)


def combine(copies: list[Copy]) -> Combined:
    """The Starr mean over one amino-acid variant's copies, and its class."""
    counts: dict[CopyState, int] = {"measured": 0, "below": 0, "above": 0, "failed": 0}
    for c in copies:
        counts[c.state] += 1
    entered = [c for c in copies if c.state != "failed"]
    n = len(entered)
    if n == 0:
        return Combined(None, None, Classification("Failed", None), None, counts)

    values = [c.log10_kd for c in entered]
    mean = math.fsum(values) / n
    sides = {c.state for c in entered}
    bound = sides == {"below"} or sides == {"above"}
    if bound:
        return Combined(mean, None, Classification("Partial", None), True, counts)

    se = _floored_se(entered) if n >= 2 else entered[0].log10_kd_se
    if se is not None:
        cls = "Good" if math.isfinite(se) and se <= KD_SD_GOOD else "Partial"
        return Combined(mean, se, Classification(cls, None), False, counts)
    return Combined(mean, None, Classification(entered[0].affinity_class, None), False, counts)


def resolve_parents(
    links: pl.DataFrame, mutations: pl.DataFrame
) -> tuple[dict[str, str], dict[str, str], list[dict[str, str]], int]:
    """Each variant's single parent, each parent's own variant, the parents that cannot be identified,
    and how many variants have several parents."""
    links = _keys(links, COL_AA, COL_PARENT)
    per_variant = links.group_by(COL_AA).agg(pl.col(COL_PARENT).unique())
    assignment: dict[str, str] = {}
    several = 0
    for aa, parents in zip(per_variant[COL_AA].to_list(), per_variant[COL_PARENT].to_list()):
        if len(parents) == 1:
            assignment[aa] = parents[0]
        else:
            several += 1

    mutation_count = pl.col(COL_MUTATIONS).cast(pl.Int64, strict=False)
    unmutated = set(mutations.filter(mutation_count == 0)[COL_AA].cast(pl.Utf8).to_list())
    candidates: dict[str, list[str]] = {}
    for aa, parent in zip(links[COL_AA].to_list(), links[COL_PARENT].to_list()):
        candidates.setdefault(parent, [])
        if aa in unmutated:
            candidates[parent].append(aa)

    own: dict[str, str] = {}
    unidentifiable: list[dict[str, str]] = []
    for parent in sorted(candidates):
        found = candidates[parent]
        if len(found) == 1:
            own[parent] = found[0]
        else:
            reason = REASON_NO_UNMUTATED if not found else REASON_SEVERAL_UNMUTATED
            unidentifiable.append({"parentId": parent, "reason": reason})
    return assignment, own, unidentifiable, several


def parent_delta(
    variant: Combined | None, parent: Combined | None, is_parent: bool
) -> tuple[float | None, float | None]:
    """deltaLog10Kd and its error. None where either side has no K_D."""
    if variant is None or parent is None or variant.log10_kd is None or parent.log10_kd is None:
        return None, None
    if is_parent:
        return 0.0, 0.0
    se = None
    if variant.log10_kd_se is not None and parent.log10_kd_se is not None:
        se = math.hypot(variant.log10_kd_se, parent.log10_kd_se)
    return variant.log10_kd - parent.log10_kd, se


def combine_synonymous(
    per_clonotype: pl.DataFrame,
    inputs: SynonymousInputs,
    *,
    log10_lo: float | None,
    log10_hi: float | None,
) -> tuple[pl.DataFrame, dict]:
    """The amino-acid frame and the run summary, from the per-nucleotide-variant fits.

    `log10_lo` and `log10_hi` are log10 of the lowest and highest non-zero concentrations: the walls
    a censored copy enters at.
    """
    fitted = per_clonotype[COL_CLONOTYPE].cast(pl.Utf8).to_list()
    mapping, unmapped, ambiguous = resolve_nt_to_aa(inputs.nt_to_aa, fitted)

    copies: dict[str, list[Copy]] = {}
    for row in per_clonotype.with_columns(pl.col(COL_CLONOTYPE).cast(pl.Utf8)).iter_rows(named=True):
        aa = mapping.get(row[COL_CLONOTYPE])
        if aa is not None:
            copies.setdefault(aa, []).append(copy_from_fit(row, log10_lo, log10_hi))
    combined = {aa: combine(copies[aa]) for aa in sorted(copies)}

    summary: dict = {
        "ntVariants": len(fitted),
        "aaVariants": len(combined),
        "unmappedNtVariants": unmapped,
        "ambiguousNtVariants": ambiguous,
        "parentMode": inputs.aa_parent_links is not None,
    }
    columns: dict[str, list] = {name: [] for name in PER_VARIANT_SCHEMA}
    for aa, c in combined.items():
        columns[COL_AA].append(aa)
        columns["kd"].append(None if c.log10_kd is None else 10.0**c.log10_kd)
        columns["kdLog10Se"].append(c.log10_kd_se)
        columns["affinityClass"].append(c.classification.affinity_class)
        columns["kdOutOfRange"].append(c.out_of_range)
        columns["ntVariants"].append(sum(c.counts.values()))
        columns["ntMeasured"].append(c.counts["measured"])
        columns["ntCensoredBelow"].append(c.counts["below"])
        columns["ntCensoredAbove"].append(c.counts["above"])
        columns["ntFailed"].append(c.counts["failed"])
    schema = dict(PER_VARIANT_SCHEMA)

    if inputs.aa_parent_links is not None and inputs.aa_mutations is not None:
        assignment, own, unidentifiable, several = resolve_parents(inputs.aa_parent_links, inputs.aa_mutations)
        deltas: list[float | None] = []
        delta_ses: list[float | None] = []
        for aa, c in combined.items():
            parent_id = assignment.get(aa)
            parent_key = own.get(parent_id) if parent_id is not None else None
            parent = combined.get(parent_key) if parent_key is not None else None
            delta, se = parent_delta(c, parent, aa == parent_key)
            deltas.append(delta)
            delta_ses.append(se)
        columns["kdLog10ParentDelta"] = deltas
        columns["kdLog10ParentDeltaSe"] = delta_ses
        schema.update(DELTA_SCHEMA)
        summary["parents"] = len(own) + len(unidentifiable)
        summary["unidentifiableParents"] = unidentifiable
        summary["multiParentVariants"] = several

    return pl.DataFrame(columns, schema=schema), summary
