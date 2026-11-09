"""Synonymous nucleotide copies combined into one value per amino-acid variant, and its parent delta.

Every nucleotide variant is fitted on its own before this module runs. Only the per-copy results are
combined here; raw counts of two copies are never pooled.

The combination is Starr's plain arithmetic mean of the copies' log10 K_D (Starr 2020, 2022):

  copies       the fitted nucleotide variants that map to exactly one amino-acid variant
  failed       a copy rated Failed, or with no usable K_D; left out of the mean and counted
  censored     a copy whose K_D sits outside the titrated range; enters at the wall it passed, and is
               counted by direction ("below" the lowest concentration, "above" the highest)
  measured     every other copy; enters at its own log10 K_D
  SE           for n >= 2 copies in the mean, the larger of two terms, the spread and the fit-error
               term: the spread is the sample SD / sqrt(n), the fit-error term sqrt(sum se_i^2) / n
               over the copies that have a fit SE. With no fit SE the SE is the spread. For n = 1, the
               copy's own fit SE where the equation reports one, else null.
  SD           the spread of the copies in the mean alone, never floored: for n >= 2, the sample SD;
               null for n = 1
  bound        no copy in the mean measured, and one side holds at least half of them and more than
               the other side; the value is that side's wall, with no SE and no SD. The replicates'
               default majority, `classify.censored_side`
  conflict     no copy in the mean measured, and neither side holds that majority: averaging
               opposite walls would land inside the range and mean nothing, so the variant fails

The class is decided once, on the combined value:

  Failed       no copy in the mean, or a conflict (`conflicting_censoring`)
  bound        the wall's class under log fluorescence; under the mean bin fit, the class the copies
               at that wall share, else Partial
  Good         SE at or below KD_SD_GOOD, with more than half the copies in the mean measured
  Partial      SE above KD_SD_GOOD or not finite, or half or more of the copies at a wall
  copy's own   n = 1 and the equation reports no SE: the combined value is that copy's fit, rated by
               its equation's own rule

Each variant gets a class, in Starr's strings. An amino-acid variant is `stop` when it is not productive,
else by its amino-acid mutation count: 0 `wildtype`, 1 `1 nonsynonymous`, more `>1 nonsynonymous`. A
nucleotide variant takes its amino-acid variant's class, except that under an unmutated amino-acid variant
it is `wildtype` with no nucleotide mutation and `synonymous` with some. Without the productive column there
is no `stop` class. A stop amino-acid variant has no K_D, SE or parent delta and is rated Failed; its
nucleotide variants keep their own fits.

The outlier test asks, per copy, whether it disagrees with the other copies of its amino-acid variant.
It runs on the final values, over an amino-acid variant whose value is a mean of k >= 2 copies:

  tested       the copies that enter that mean, measured and censored alike (a censored copy at its
               wall); a failed copy is not. A bound, a failed and a stop amino-acid variant have no
               mean, so none of its copies is tested
  delta        the copy's log10 K_D minus the mean of the other copies'
  sigma        sqrt(se^2 + se_others^2), where the copy has a fit SE and the others' SE exists: the
               others' SE by the SE rule above over n = k - 1, the other copy's own fit SE for n = 1.
               Otherwise the dataset sigma * sqrt(1 + 1/(k - 1))
  dataset      1.4826 * MAD of every tested copy's deviation from its own amino-acid mean, pooled
  sigma        over the run; robust to the few copies that are real outliers
  z            delta / sigma; flagged at |z| >= SYNONYM_OUTLIER_Z. Null where sigma is null or 0

With k = 2 the test is symmetric: both copies get the same |z|, and the amino-acid note names a
disagreement, not an odd copy.

The parent delta is log10 K_D of the variant minus that of its parent, both at amino-acid grain, with
error sqrt(se_v^2 + se_p^2). A parent's own row is 0 with error 0. A parent is the one amino-acid variant
with zero mutations linked to it; a parent with none or several cannot be identified. The nucleotide grain
has its own parent delta, by the same rule on the nucleotide parent link and mutation count: each fitted
nt variant's own log10 K_D, as it enters the mean, minus that of the one unmutated nt variant of its parent.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Literal, Protocol

import polars as pl

from classify import REASON_CONFLICTING, Classification, bound_class, censored_side
from constants import COL_CLONOTYPE, DEFAULT_CENSORED_MAJORITY, KD_SD_GOOD, MAD_TO_SD, SYNONYM_OUTLIER_Z

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
# The copies behind one replicate's amino-acid value, carried into the per-replicate rows: the
# spread of the measured copies inside that replicate, and how many copies there were.
PER_REPLICATE_COPY_COLUMNS: dict[str, pl.DataType] = {
    "kdLog10SynonymSd": pl.Float64,
    "ntVariants": pl.Int64,
    "ntMeasured": pl.Int64,
}
CLASS_SCHEMA: dict[str, pl.DataType] = {COL_CLASS: pl.Utf8}
DELTA_SCHEMA: dict[str, pl.DataType] = {
    "kdLog10ParentDelta": pl.Float64,
    "kdLog10ParentDeltaSe": pl.Float64,
}
# The outlier test, per nt variant; null for a copy that is not tested.
OUTLIER_SCHEMA: dict[str, pl.DataType] = {
    "kdLog10SynonymDelta": pl.Float64,
    "kdLog10SynonymDeltaSe": pl.Float64,
    "kdLog10SynonymZ": pl.Float64,
    "synonymOutlier": pl.Boolean,
}
# Per amino-acid variant: the flagged copies, in words; null when none is flagged.
COL_OUTLIERS = "synonymOutliers"


@dataclass(frozen=True)
class SynonymousInputs:
    """The nt -> aa linker and the per-variant columns that resolved.

    The amino-acid parent delta needs `aa_mutations` and `aa_parent_links`; the amino-acid class needs
    `aa_mutations`; the nucleotide class needs `aa_mutations` and `nt_mutations`; `stop` needs
    `aa_productive`.
    """

    nt_to_aa: pl.DataFrame  # [ntKey, aaKey]
    aa_mutations: pl.DataFrame | None = None  # [aaKey, mutationCount]
    aa_parent_links: pl.DataFrame | None = None  # [aaKey, parentId]
    aa_productive: pl.DataFrame | None = None  # [aaKey, productive], "true" / "false"
    nt_mutations: pl.DataFrame | None = None  # [ntKey, mutationCount]
    nt_parent_links: pl.DataFrame | None = None  # [ntKey, parentId]

    def __post_init__(self) -> None:
        if self.aa_parent_links is not None and self.aa_mutations is None:
            raise ValueError("the amino-acid parent links need the amino-acid mutation counts")
        if self.nt_parent_links is not None and self.nt_mutations is None:
            raise ValueError("the nucleotide parent links need the nucleotide mutation counts")


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


COL_LABEL = "label"
COL_AA_VARIANT = "aaVariant"


def nt_aa_variants(
    nt_to_aa: pl.DataFrame, aa_labels: pl.DataFrame | None, curved: Collection[str] | None = None
) -> pl.DataFrame:
    """Each nt variant's amino-acid variant: its key, and its label where the graphs tab by it.

    One row per nt variant that maps to exactly one amino-acid variant, with that variant's key.
    The label (else the key) is set only on an nt variant in `curved`, and only where its
    amino-acid variant has at least 2 of them: the synonymous graphs list exactly those, on the
    nt axis, where the linker cannot reach. `curved` None counts every nt variant.
    """
    per_nt = _keys(nt_to_aa, COL_NT, COL_AA).group_by(COL_NT).agg(pl.col(COL_AA)).filter(pl.col(COL_AA).list.len() == 1)
    frame = per_nt.select(pl.col(COL_NT).alias(COL_CLONOTYPE), pl.col(COL_AA).list.first())
    has_curve = pl.lit(True) if curved is None else pl.col(COL_CLONOTYPE).is_in(list(curved))
    frame = frame.with_columns(has_curve.alias("_curved"))
    frame = frame.with_columns((pl.col("_curved") & (pl.col("_curved").sum().over(COL_AA) >= 2)).alias("_shown"))
    if aa_labels is not None:
        labels = _keys(aa_labels, COL_AA, COL_LABEL).unique(subset=[COL_AA], keep="first")
        frame = frame.join(labels, on=COL_AA, how="left")
    else:
        frame = frame.with_columns(pl.lit(None, dtype=pl.Utf8).alias(COL_LABEL))
    variant = pl.when(pl.col("_shown")).then(pl.coalesce(COL_LABEL, COL_AA)).otherwise(None)
    return frame.select(COL_CLONOTYPE, COL_AA, variant.alias(COL_AA_VARIANT)).sort(COL_CLONOTYPE)


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


def _floored_se(entered: Sequence[HasKd]) -> float:
    """max(spread, fit-error term) over n >= 2 copies; the spread alone when no copy has a
    fit SE."""
    n = len(entered)
    spread = statistics.stdev([c.log10_kd for c in entered]) / math.sqrt(n)
    fit_ses = [c.log10_kd_se for c in entered if c.log10_kd_se is not None]
    if not fit_ses:
        return spread
    return max(spread, math.sqrt(math.fsum(se * se for se in fit_ses)) / n)


def combine(copies: list[Copy], *, wall_classes: bool = True) -> Combined:
    """The mean over one amino-acid variant's copies, and its class. `wall_classes` is True under the
    equation that rates a bound by its wall."""
    counts: dict[CopyState, int] = {"measured": 0, "below": 0, "above": 0, "failed": 0}
    for c in copies:
        counts[c.state] += 1
    entered = [c for c in copies if c.state != "failed"]
    n = len(entered)
    if n == 0:
        return Combined(None, None, Classification("Failed", None), None, counts)

    mean = math.fsum(c.log10_kd for c in entered) / n
    if not any(c.state == "measured" for c in entered):
        side = censored_side([c.state for c in entered], DEFAULT_CENSORED_MAJORITY)
        if side is None:
            return Combined(None, None, Classification("Failed", REASON_CONFLICTING), None, counts)
        on_side = [c for c in entered if c.state == side]
        cls = bound_class(side, [c.affinity_class for c in on_side], wall_classes=wall_classes)
        return Combined(on_side[0].log10_kd, None, Classification(cls, None), True, counts)

    if n == 1:
        se, sd = entered[0].log10_kd_se, None
    else:
        se = _floored_se(entered)
        sd = statistics.stdev([c.log10_kd for c in entered])
    if se is not None:
        # A copy at a wall enters at the wall, so two of them agree exactly however wide the truth
        # is. Good needs more than half the copies in the mean measured.
        measured_majority = 2 * sum(c.state == "measured" for c in entered) > n
        cls = "Good" if math.isfinite(se) and se <= KD_SD_GOOD and measured_majority else "Partial"
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


def _is_bound(value: HasKd) -> bool:
    """True for a value censored at a wall of the titration."""
    return bool(getattr(value, "out_of_range", False)) or getattr(value, "state", None) in ("below", "above")


def synonym_spread(
    aa_keys: list[str], nt_values: Mapping[str, HasKd], mapping: dict[str, str], stops: set[str]
) -> list[float | None]:
    """Per amino-acid variant, in `aa_keys` order: the SD of log10 Kd,app over its synonymous nt
    variants measured inside the range. None with fewer than two, and for a stop variant."""
    measured: dict[str, list[float]] = {}
    for nt, value in nt_values.items():
        aa = mapping.get(nt)
        if aa is None or aa in stops or value.log10_kd is None or _is_bound(value):
            continue
        measured.setdefault(aa, []).append(value.log10_kd)
    return [statistics.stdev(measured[aa]) if len(measured.get(aa, [])) >= 2 else None for aa in aa_keys]


@dataclass(frozen=True)
class OutlierTest:
    """One copy against the other copies of its amino-acid variant. For a copy at a wall (`bound`),
    `delta` and `z` are the least the gap can be."""

    delta: float
    sigma: float | None
    z: float | None
    bound: bool = False

    @property
    def flagged(self) -> bool | None:
        return None if self.z is None else abs(self.z) >= SYNONYM_OUTLIER_Z


def outlier_groups(
    aa_values: Mapping[str, HasKd], nt_values: Mapping[str, HasKd], mapping: dict[str, str], stops: set[str]
) -> dict[str, list[tuple[str, HasKd]]]:
    """Per amino-acid variant whose value is a mean of at least two copies: the copies in that mean,
    by nt key. A copy enters when it has a log10 K_D, its own or its wall's."""
    entered: dict[str, list[tuple[str, HasKd]]] = {}
    for nt in sorted(nt_values):
        aa = mapping.get(nt)
        if aa is not None and nt_values[nt].log10_kd is not None:
            entered.setdefault(aa, []).append((nt, nt_values[nt]))
    groups: dict[str, list[tuple[str, HasKd]]] = {}
    for aa, value in aa_values.items():
        if aa in stops or value.log10_kd is None or _is_bound(value) or len(entered.get(aa, [])) < 2:
            continue
        groups[aa] = entered[aa]
    return groups


def _deviations(groups: Mapping[str, list[tuple[str, HasKd]]]) -> list[float]:
    """Each copy measured inside the range minus the mean of the m >= 2 measured copies of its
    amino-acid variant, times sqrt(m / (m - 1)): the mean includes the copy, which shrinks the
    deviation's SD to sigma * sqrt((m - 1) / m). A copy at a wall is left out: its value is a bound."""
    out: list[float] = []
    for copies in groups.values():
        values = [c.log10_kd for _, c in copies if not _is_bound(c)]
        m = len(values)
        if m < 2:
            continue
        mean = math.fsum(values) / m
        inflate = math.sqrt(m / (m - 1))
        out.extend((v - mean) * inflate for v in values)
    return out


def dataset_sigma(groups: Mapping[str, list[tuple[str, HasKd]]]) -> float | None:
    """1.4826 * MAD of the pooled deviations; None with no deviation."""
    deviations = _deviations(groups)
    if not deviations:
        return None
    center = statistics.median(deviations)
    return MAD_TO_SD * statistics.median(abs(d - center) for d in deviations)


def outlier_tests(copies: list[tuple[str, HasKd]], sigma: float | None) -> dict[str, OutlierTest]:
    """Each copy of one amino-acid variant against the mean of its other copies measured inside the
    range. A copy at a wall is never among the others: its value is a bound, and would drag their
    mean. A copy with no other measured copy is not tested. `sigma` is the dataset sigma, used where
    a copy or the others have no SE."""
    measured = [(nt, c) for nt, c in copies if not _is_bound(c)]
    tests: dict[str, OutlierTest] = {}
    for nt, copy in copies:
        others = [c for other, c in measured if other != nt]
        m = len(others)
        if m == 0:
            continue
        delta = copy.log10_kd - math.fsum(c.log10_kd for c in others) / m
        se_others = others[0].log10_kd_se if m == 1 else _floored_se(others)
        if copy.log10_kd_se is not None and se_others is not None:
            scale = math.hypot(copy.log10_kd_se, se_others)
        elif sigma is not None:
            scale = sigma * math.sqrt(1 + 1 / m)
        else:
            scale = None
        z = None if scale is None or not math.isfinite(scale) or scale <= 0 else delta / scale
        tests[nt] = OutlierTest(delta, scale, z, _is_bound(copy))
    return tests


def outlier_note(tests: Mapping[str, OutlierTest]) -> str | None:
    """The flagged copies of one amino-acid variant, in words; None when none is flagged.

    With two copies measured inside the range both are flagged or neither, so the note names the
    pair. Otherwise each flagged copy, largest |z| first: weaker for a higher log10 K_D, tighter for
    a lower one. A copy at a wall reads as the least gap, without a z.
    """
    flagged = sorted(((nt, t) for nt, t in tests.items() if t.flagged), key=lambda p: (-abs(p[1].z), p[0]))
    if not flagged:
        return None
    if len(tests) == 2 and not any(t.bound for t in tests.values()):
        a, b = sorted(tests)
        t = tests[a]
        return f"{a} vs {b}: disagree, {abs(t.delta):.2f} decades apart, |z| {abs(t.z):.2f}"

    def describe(nt: str, t: OutlierTest) -> str:
        side = "weaker" if t.delta > 0 else "tighter"
        if t.bound:
            return f"{nt}: {side}, at least {abs(t.delta):.2f} decades, {'above' if t.delta > 0 else 'below'} range"
        return f"{nt}: {side}, {t.delta:+.2f} decades, z {t.z:+.2f}"

    return "; ".join(describe(nt, t) for nt, t in flagged)


def parent_deltas(
    values: Mapping[str, HasKd], links: pl.DataFrame, mutations: pl.DataFrame, key: str
) -> tuple[list[float | None], list[float | None], dict]:
    """The parent delta and its error for each variant of `values`, in its order, and the summary
    counts."""
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
    wall_classes: bool = True,
) -> tuple[pl.DataFrame, pl.DataFrame, dict]:
    """The amino-acid frame, the nucleotide frame, and the run summary, from the per-nt-variant fits.

    The nucleotide frame holds one row per fitted nt variant: the outlier test, and the nt-grain class
    and parent delta where their inputs came.

    `log10_lo` and `log10_hi` are log10 of the lowest and highest non-zero concentrations: the walls
    a censored copy enters at.
    """
    fitted = per_clonotype[COL_CLONOTYPE].cast(pl.Utf8).to_list()
    mapping, unmapped, ambiguous = resolve_nt_to_aa(inputs.nt_to_aa, fitted)

    rows = list(per_clonotype.with_columns(pl.col(COL_CLONOTYPE).cast(pl.Utf8)).iter_rows(named=True))
    copies: dict[str, list[Copy]] = {}
    for row in rows:
        aa = mapping.get(row[COL_CLONOTYPE])
        if aa is not None:
            copies.setdefault(aa, []).append(copy_from_fit(row, log10_lo, log10_hi))
    combined = {aa: combine(copies[aa], wall_classes=wall_classes) for aa in sorted(copies)}
    stops = stop_variants(inputs.aa_productive) if inputs.aa_productive is not None else set()
    # A stop variant has no value: its counts stay, every value goes.
    combined = {aa: replace(STOP, counts=c.counts) if aa in stops else c for aa, c in combined.items()}

    summary: dict = {
        "ntVariants": len(fitted),
        "aaVariants": len(combined),
        "unmappedNtVariants": unmapped,
        "ambiguousNtVariants": ambiguous,
        "parentMode": inputs.aa_parent_links is not None,
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
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The class, parent-delta and outlier columns of the amino-acid frame, and the nucleotide frame.

    Decided once, on the final values: `aa_values` in the order of `per_variant`'s rows, `nt_values`
    one per fitted nt variant. Adds the parent and outlier counts and the dataset sigma to `summary`.
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

    spread = synonym_spread(list(aa_values), nt_values, mapping, stops)
    per_variant = per_variant.with_columns(pl.Series("kdLog10SynonymSd", spread, dtype=pl.Float64))
    present = [v for v in spread if v is not None]
    summary["synonymSpread"] = {"variants": len(present), "median": statistics.median(present) if present else None}

    groups = outlier_groups(aa_values, nt_values, mapping, stops)
    sigma = dataset_sigma(groups)
    per_aa_tests = {aa: outlier_tests(copies, sigma) for aa, copies in groups.items()}
    tests = {nt: t for group in per_aa_tests.values() for nt, t in group.items()}
    notes = [outlier_note(per_aa_tests[aa]) if aa in per_aa_tests else None for aa in aa_values]
    per_variant = per_variant.with_columns(pl.Series(COL_OUTLIERS, notes, dtype=pl.Utf8))
    summary["synonymOutlier"] = {
        "zThreshold": SYNONYM_OUTLIER_Z,
        "datasetSigma": sigma,
        "testedAaVariants": len(groups),
        "testedNtVariants": len(tests),
        "flaggedNtVariants": sum(1 for t in tests.values() if t.flagged),
        "flaggedAaVariants": sum(1 for n in notes if n is not None),
    }

    fitted = list(nt_values)
    nt_columns: dict[str, list] = {COL_NT: fitted}
    nt_schema: dict[str, pl.DataType] = {COL_NT: pl.Utf8}
    if aa_mutations is not None and inputs.nt_mutations is not None:
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
    nt_tests = [tests.get(nt) for nt in fitted]
    nt_columns["kdLog10SynonymDelta"] = [None if t is None else t.delta for t in nt_tests]
    nt_columns["kdLog10SynonymDeltaSe"] = [None if t is None else t.sigma for t in nt_tests]
    nt_columns["kdLog10SynonymZ"] = [None if t is None else t.z for t in nt_tests]
    nt_columns["synonymOutlier"] = [None if t is None else t.flagged for t in nt_tests]
    nt_schema.update(OUTLIER_SCHEMA)
    return per_variant, pl.DataFrame(nt_columns, schema=nt_schema)
