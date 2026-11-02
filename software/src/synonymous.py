"""Synonymous nucleotide copies combined into one value per amino-acid variant, and its parent delta.

Every nucleotide variant is fitted on its own before this module runs. Only the per-copy results are
combined here; raw counts of two copies are never pooled.

The combination is the mean of the copies' log10 K_D. Unweighted, it is Starr's plain arithmetic mean
(Starr 2020, 2022). Read-weighted, each copy carries its total reads, w_i, over the samples that enter the
fit; a copy with no reads carries no weight and stays out of the mean:

  copies       the fitted nucleotide variants that map to exactly one amino-acid variant
  failed       a copy rated Failed, or with no usable K_D; left out of the mean and counted
  censored     a copy whose K_D sits outside the titrated range; enters at the wall it passed, and is
               counted by direction ("below" the lowest concentration, "above" the highest)
  measured     every other copy; enters at its own log10 K_D
  SE           for n >= 2 copies in the mean, the larger of two terms, the spread and the fit-error
               term. The fit-error term sums over the copies that have a fit SE; it is null when none
               has one, and the SE is then the spread. For n = 1, the copy's own fit SE where the
               equation reports one, else null.
                 unweighted   spread sample SD / sqrt(n); fit-error sqrt(sum se_i^2) / n
                 weighted     spread sqrt(n / (n - 1) * sum w_i^2 (x_i - x_w)^2) / sum w_i, around the
                              weighted mean x_w; fit-error sqrt(sum w_i^2 se_i^2) / sum w_i
  SD           the spread of the copies in the mean alone, never floored: for n >= 2, the sample SD
               unweighted, sqrt(n / (n - 1) * sum w_i (x_i - x_w)^2 / sum w_i) weighted; null for n = 1
  bound        every copy in the mean censored on the same side; the value is that wall, with no SE
               and no SD

The class is decided once, on the combined value:

  Failed       no copy in the mean
  Partial      a bound
  Good         SE at or below KD_SD_GOOD
  Partial      SE above KD_SD_GOOD, or not finite
  copy's own   n = 1 and the equation reports no SE: the combined value is that copy's fit, rated by
               its equation's own rule

Each variant gets a class, in Starr's strings. An amino-acid variant is `stop` when it is not productive,
else by its amino-acid mutation count: 0 `wildtype`, 1 `1 nonsynonymous`, more `>1 nonsynonymous`. A
nucleotide variant takes its amino-acid variant's class, except that under an unmutated amino-acid variant
it is `wildtype` with no nucleotide mutation and `synonymous` with some. Without the productive column there
is no `stop` class. A stop amino-acid variant has no K_D, SE or parent delta and is rated Failed; its
nucleotide variants keep their own fits.

The parent delta is log10 K_D of the variant minus that of its parent, both at amino-acid grain, with
error sqrt(se_v^2 + se_p^2). A parent's own row is 0 with error 0. A parent is the one amino-acid variant
with zero mutations linked to it; a parent with none or several cannot be identified. The nucleotide grain
has its own parent delta, by the same rule on the nucleotide parent link and mutation count: each fitted
nt variant's own log10 K_D, as it enters the mean, minus that of the one unmutated nt variant of its parent.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import Literal, Protocol

import polars as pl

from classify import Classification, bound_class
from constants import COL_CLONOTYPE, KD_SD_GOOD

COL_NT = "ntKey"
COL_AA = "aaKey"
COL_PARENT = "parentId"
COL_MUTATIONS = "mutationCount"

COL_PRODUCTIVE = "productive"
COL_CLASS = "variantClass"

CLASS_WILDTYPE = "wildtype"
CLASS_SYNONYMOUS = "synonymous"
CLASS_STOP = "stop"
CLASS_ONE = "1 nonsynonymous"
CLASS_SEVERAL = ">1 nonsynonymous"

REASON_NO_UNMUTATED = "no-unmutated-variant"
REASON_SEVERAL_UNMUTATED = "several-unmutated-variants"

CopyState = Literal["measured", "below", "above", "failed"]
Weighting = Literal["none", "reads"]
WEIGHTINGS: tuple[Weighting, ...] = ("none", "reads")

PER_VARIANT_SCHEMA: dict[str, pl.DataType] = {
    COL_AA: pl.Utf8,
    "kd": pl.Float64,
    "kdNegLog10": pl.Float64,
    "kdLog10Se": pl.Float64,
    "kdLog10Sd": pl.Float64,
    "affinityClass": pl.Utf8,
    "kdOutOfRange": pl.Boolean,
    "ntVariants": pl.Int64,
    "ntMeasured": pl.Int64,
    "ntCensoredBelow": pl.Int64,
    "ntCensoredAbove": pl.Int64,
    "ntFailed": pl.Int64,
}
CLASS_SCHEMA: dict[str, pl.DataType] = {COL_CLASS: pl.Utf8}
DELTA_SCHEMA: dict[str, pl.DataType] = {
    "kdLog10ParentDelta": pl.Float64,
    "kdLog10ParentDeltaSe": pl.Float64,
}


@dataclass(frozen=True)
class SynonymousInputs:
    """The nt -> aa linker, the per-variant columns that resolved, and the weighting.

    The amino-acid parent delta needs `aa_mutations` and `aa_parent_links`; the amino-acid class needs
    `aa_mutations`; the nucleotide class needs `aa_mutations` and `nt_mutations`; `stop` needs
    `aa_productive`.
    """

    nt_to_aa: pl.DataFrame  # [ntKey, aaKey]
    aa_mutations: pl.DataFrame | None = None  # [aaKey, mutationCount]
    aa_parent_links: pl.DataFrame | None = None  # [aaKey, parentId]
    weighting: Weighting = "none"
    aa_productive: pl.DataFrame | None = None  # [aaKey, productive], "true" / "false"
    nt_mutations: pl.DataFrame | None = None  # [ntKey, mutationCount]
    nt_parent_links: pl.DataFrame | None = None  # [ntKey, parentId]

    def __post_init__(self) -> None:
        if self.aa_parent_links is not None and self.aa_mutations is None:
            raise ValueError("the amino-acid parent links need the amino-acid mutation counts")
        if self.nt_parent_links is not None and self.nt_mutations is None:
            raise ValueError("the nucleotide parent links need the nucleotide mutation counts")
        if self.weighting not in WEIGHTINGS:
            raise ValueError(f"unknown weighting: {self.weighting}")


class HasKd(Protocol):
    """A value the parent delta reads: a copy, a synonymous combination, or a replicate mean."""

    @property
    def log10_kd(self) -> float | None: ...

    @property
    def log10_kd_se(self) -> float | None: ...


@dataclass(frozen=True)
class Copy:
    state: CopyState
    log10_kd: float | None
    log10_kd_se: float | None
    affinity_class: str
    weight: float = 1.0


@dataclass(frozen=True)
class Combined:
    log10_kd: float | None
    log10_kd_se: float | None
    classification: Classification
    out_of_range: bool | None
    counts: dict[CopyState, int]
    log10_kd_sd: float | None = None


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
    """max(spread, fit-error term) over n >= 2 unweighted copies; the spread alone when no copy has a
    fit SE."""
    n = len(entered)
    spread = statistics.stdev([c.log10_kd for c in entered]) / math.sqrt(n)
    fit_ses = [c.log10_kd_se for c in entered if c.log10_kd_se is not None]
    if not fit_ses:
        return spread
    return max(spread, math.sqrt(math.fsum(se * se for se in fit_ses)) / n)


def _floored_weighted_se(entered: list[Copy], mean: float) -> float:
    """The read-weighted counterpart of `_floored_se`, around the weighted mean."""
    n = len(entered)
    total = math.fsum(c.weight for c in entered)
    scatter = math.fsum(c.weight**2 * (c.log10_kd - mean) ** 2 for c in entered)
    spread = math.sqrt(n / (n - 1) * scatter) / total
    fit = [c.weight**2 * c.log10_kd_se**2 for c in entered if c.log10_kd_se is not None]
    if not fit:
        return spread
    return max(spread, math.sqrt(math.fsum(fit)) / total)


def _spread_sd(entered: list[Copy], mean: float, weighted: bool) -> float:
    """The SD of n >= 2 copies around `mean`."""
    if not weighted:
        return statistics.stdev([c.log10_kd for c in entered])
    n = len(entered)
    total = math.fsum(c.weight for c in entered)
    scatter = math.fsum(c.weight * (c.log10_kd - mean) ** 2 for c in entered)
    return math.sqrt(n / (n - 1) * scatter / total)


def combine(copies: list[Copy], *, weighted: bool = False, wall_classes: bool = True) -> Combined:
    """The mean over one amino-acid variant's copies, and its class. `weighted` reads each copy's
    `weight`. `wall_classes` is True under the equation that rates a bound by its wall."""
    counts: dict[CopyState, int] = {"measured": 0, "below": 0, "above": 0, "failed": 0}
    for c in copies:
        counts[c.state] += 1
    entered = [c for c in copies if c.state != "failed" and (not weighted or c.weight > 0)]
    n = len(entered)
    if n == 0:
        return Combined(None, None, Classification("Failed", None), None, counts)

    if weighted:
        mean = math.fsum(c.weight * c.log10_kd for c in entered) / math.fsum(c.weight for c in entered)
    else:
        mean = math.fsum(c.log10_kd for c in entered) / n
    sides = {c.state for c in entered}
    bound = sides == {"below"} or sides == {"above"}
    if bound:
        side = entered[0].state
        cls = bound_class(side, [c.affinity_class for c in entered], wall_classes=wall_classes)
        return Combined(mean, None, Classification(cls, None), True, counts)

    if n == 1:
        se, sd = entered[0].log10_kd_se, None
    else:
        se = _floored_weighted_se(entered, mean) if weighted else _floored_se(entered)
        sd = _spread_sd(entered, mean, weighted)
    if se is not None:
        cls = "Good" if math.isfinite(se) and se <= KD_SD_GOOD else "Partial"
        return Combined(mean, se, Classification(cls, None), False, counts, sd)
    return Combined(mean, None, Classification(entered[0].affinity_class, None), False, counts, sd)


def mutation_counts(frame: pl.DataFrame, key: str) -> dict[str, int]:
    counts = frame.select(pl.col(key).cast(pl.Utf8), pl.col(COL_MUTATIONS).cast(pl.Int64, strict=False)).drop_nulls()
    return dict(zip(counts[key].to_list(), counts[COL_MUTATIONS].to_list()))


def stop_variants(productive: pl.DataFrame) -> set[str]:
    """The amino-acid variants marked not productive."""
    flag = pl.col(COL_PRODUCTIVE).cast(pl.Utf8).str.strip_chars().str.to_lowercase()
    return set(productive.filter(flag == "false")[COL_AA].cast(pl.Utf8).to_list())


def aa_class(mutations: int | None, stop: bool) -> str | None:
    """Starr's class of an amino-acid variant; None when its mutation count is unknown."""
    if stop:
        return CLASS_STOP
    if mutations is None:
        return None
    if mutations == 0:
        return CLASS_WILDTYPE
    return CLASS_ONE if mutations == 1 else CLASS_SEVERAL


def nt_class(aa_mutations: int | None, nt_mutations: int | None, stop: bool) -> str | None:
    """Starr's class of a nucleotide variant, from its amino-acid variant; None when unknown."""
    cls = aa_class(aa_mutations, stop)
    if cls != CLASS_WILDTYPE:
        return cls
    if nt_mutations is None:
        return None
    return CLASS_WILDTYPE if nt_mutations == 0 else CLASS_SYNONYMOUS


STOP = Combined(None, None, Classification("Failed", None), None, {})


def resolve_parents(
    links: pl.DataFrame, mutations: pl.DataFrame, key: str = COL_AA
) -> tuple[dict[str, str], dict[str, str], list[dict[str, str]], int]:
    """Each variant's single parent, each parent's own variant, the parents that cannot be identified,
    and how many variants have several parents. `key` names the variant column of both frames."""
    links = _keys(links, key, COL_PARENT)
    per_variant = links.group_by(key).agg(pl.col(COL_PARENT).unique())
    assignment: dict[str, str] = {}
    several = 0
    for variant, parents in zip(per_variant[key].to_list(), per_variant[COL_PARENT].to_list()):
        if len(parents) == 1:
            assignment[variant] = parents[0]
        else:
            several += 1

    mutation_count = pl.col(COL_MUTATIONS).cast(pl.Int64, strict=False)
    unmutated = set(mutations.filter(mutation_count == 0)[key].cast(pl.Utf8).to_list())
    candidates: dict[str, list[str]] = {}
    for variant, parent in zip(links[key].to_list(), links[COL_PARENT].to_list()):
        candidates.setdefault(parent, [])
        if variant in unmutated:
            candidates[parent].append(variant)

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


def parent_delta(variant: HasKd | None, parent: HasKd | None, is_parent: bool) -> tuple[float | None, float | None]:
    """deltaLog10Kd and its error. None where either side has no K_D."""
    if variant is None or parent is None or variant.log10_kd is None or parent.log10_kd is None:
        return None, None
    if is_parent:
        return 0.0, 0.0
    se = None
    if variant.log10_kd_se is not None and parent.log10_kd_se is not None:
        se = math.hypot(variant.log10_kd_se, parent.log10_kd_se)
    return variant.log10_kd - parent.log10_kd, se


def parent_deltas(
    values: Mapping[str, HasKd], links: pl.DataFrame, mutations: pl.DataFrame, key: str
) -> tuple[list[float | None], list[float | None], dict]:
    """The parent delta and its error for each variant of `values`, in its order, and the summary counts."""
    assignment, own, unidentifiable, several = resolve_parents(links, mutations, key)
    deltas: list[float | None] = []
    ses: list[float | None] = []
    for variant, value in values.items():
        parent_id = assignment.get(variant)
        parent_key = own.get(parent_id) if parent_id is not None else None
        parent = values.get(parent_key) if parent_key is not None else None
        delta, se = parent_delta(value, parent, variant == parent_key)
        deltas.append(delta)
        ses.append(se)
    counts = {
        "parents": len(own) + len(unidentifiable),
        "unidentifiableParents": unidentifiable,
        "multiParentVariants": several,
    }
    return deltas, ses, counts


def combine_synonymous(
    per_clonotype: pl.DataFrame,
    inputs: SynonymousInputs,
    *,
    log10_lo: float | None,
    log10_hi: float | None,
    reads: dict[str, float] | None = None,
    wall_classes: bool = True,
) -> tuple[pl.DataFrame, pl.DataFrame | None, dict]:
    """The amino-acid frame, the nucleotide frame, and the run summary, from the per-nt-variant fits.

    The nucleotide frame, one row per fitted nt variant, exists when `nt_mutations` comes; it holds
    the nt-grain columns that resolved.

    `log10_lo` and `log10_hi` are log10 of the lowest and highest non-zero concentrations: the walls
    a censored copy enters at. `reads` is each nt variant's total reads over the samples that enter
    the fit; read-weighting needs it, and a variant absent from it carries no weight.
    """
    weighted = inputs.weighting == "reads"
    if weighted and reads is None:
        raise ValueError("read-weighting needs the reads of each nt variant")
    fitted = per_clonotype[COL_CLONOTYPE].cast(pl.Utf8).to_list()
    mapping, unmapped, ambiguous = resolve_nt_to_aa(inputs.nt_to_aa, fitted)

    rows = list(per_clonotype.with_columns(pl.col(COL_CLONOTYPE).cast(pl.Utf8)).iter_rows(named=True))
    copies: dict[str, list[Copy]] = {}
    for row in rows:
        aa = mapping.get(row[COL_CLONOTYPE])
        if aa is not None:
            c = copy_from_fit(row, log10_lo, log10_hi)
            if weighted:
                c = replace(c, weight=float(reads.get(row[COL_CLONOTYPE], 0.0)))
            copies.setdefault(aa, []).append(c)
    combined = {aa: combine(copies[aa], weighted=weighted, wall_classes=wall_classes) for aa in sorted(copies)}
    stops = stop_variants(inputs.aa_productive) if inputs.aa_productive is not None else set()
    # A stop variant has no value: its counts stay, every value goes.
    combined = {aa: replace(STOP, counts=c.counts) if aa in stops else c for aa, c in combined.items()}

    summary: dict = {
        "ntVariants": len(fitted),
        "aaVariants": len(combined),
        "unmappedNtVariants": unmapped,
        "ambiguousNtVariants": ambiguous,
        "parentMode": inputs.aa_parent_links is not None,
        "weighting": inputs.weighting,
    }
    if inputs.aa_productive is not None:
        summary["stopAaVariants"] = sum(1 for aa in combined if aa in stops)
    columns: dict[str, list] = {name: [] for name in PER_VARIANT_SCHEMA}
    for aa, c in combined.items():
        columns[COL_AA].append(aa)
        columns["kd"].append(None if c.log10_kd is None else 10.0**c.log10_kd)
        columns["kdNegLog10"].append(None if c.log10_kd is None else -c.log10_kd)
        columns["kdLog10Se"].append(c.log10_kd_se)
        columns["kdLog10Sd"].append(c.log10_kd_sd)
        columns["affinityClass"].append(c.classification.affinity_class)
        columns["kdOutOfRange"].append(c.out_of_range)
        columns["ntVariants"].append(sum(c.counts.values()))
        columns["ntMeasured"].append(c.counts["measured"])
        columns["ntCensoredBelow"].append(c.counts["below"])
        columns["ntCensoredAbove"].append(c.counts["above"])
        columns["ntFailed"].append(c.counts["failed"])
    nt_copies = {row[COL_CLONOTYPE]: copy_from_fit(row, log10_lo, log10_hi) for row in rows}
    per_variant, per_nt = annotate_variants(
        pl.DataFrame(columns, schema=PER_VARIANT_SCHEMA), combined, nt_copies, mapping, stops, inputs, summary
    )
    return per_variant, per_nt, summary


def annotate_variants(
    per_variant: pl.DataFrame,
    aa_values: Mapping[str, HasKd],
    nt_values: Mapping[str, HasKd],
    mapping: dict[str, str],
    stops: set[str],
    inputs: SynonymousInputs,
    summary: dict,
) -> tuple[pl.DataFrame, pl.DataFrame | None]:
    """The class and parent-delta columns of the amino-acid frame, and the nucleotide frame.

    Decided once, on the final values: `aa_values` in the order of `per_variant`'s rows, `nt_values`
    one per fitted nt variant. Adds the parent counts to `summary`.
    """
    aa_mutations = mutation_counts(inputs.aa_mutations, COL_AA) if inputs.aa_mutations is not None else None
    if aa_mutations is not None:
        per_variant = per_variant.with_columns(
            pl.Series(COL_CLASS, [aa_class(aa_mutations.get(aa), aa in stops) for aa in aa_values], dtype=pl.Utf8)
        )

    if inputs.aa_parent_links is not None and inputs.aa_mutations is not None:
        deltas, delta_ses, counts = parent_deltas(aa_values, inputs.aa_parent_links, inputs.aa_mutations, COL_AA)
        per_variant = per_variant.with_columns(
            pl.Series("kdLog10ParentDelta", deltas, dtype=pl.Float64),
            pl.Series("kdLog10ParentDeltaSe", delta_ses, dtype=pl.Float64),
        )
        summary.update(counts)

    if inputs.nt_mutations is None:
        return per_variant, None
    fitted = list(nt_values)
    nt_columns: dict[str, list] = {COL_NT: fitted}
    nt_schema: dict[str, pl.DataType] = {COL_NT: pl.Utf8}
    if aa_mutations is not None:
        nt_mutations = mutation_counts(inputs.nt_mutations, COL_NT)
        nt_columns[COL_CLASS] = [
            None
            if nt not in mapping
            else nt_class(aa_mutations.get(mapping[nt]), nt_mutations.get(nt), mapping[nt] in stops)
            for nt in fitted
        ]
        nt_schema.update(CLASS_SCHEMA)
    summary["ntParentMode"] = inputs.nt_parent_links is not None
    if inputs.nt_parent_links is not None:
        deltas, delta_ses, counts = parent_deltas(nt_values, inputs.nt_parent_links, inputs.nt_mutations, COL_NT)
        nt_columns["kdLog10ParentDelta"] = deltas
        nt_columns["kdLog10ParentDeltaSe"] = delta_ses
        nt_schema.update(DELTA_SCHEMA)
        summary.update({"nt" + k[0].upper() + k[1:]: v for k, v in counts.items()})
    return per_variant, pl.DataFrame(nt_columns, schema=nt_schema)
