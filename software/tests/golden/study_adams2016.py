"""Adams 2016 as the block sees it: the authors' read counts per nucleotide variant, three replicates,
four gates with their sort fractions and measured log10 fluorescence, run through `pipeline.run` with
the block's default settings.

The gates were set anew on each sort day, so each replicate's gates sit at their own fluorescence:
gate 3 at 3.35, 3.10 and 3.17 log10 in replicates 1, 2 and 3. The authors' fit reads each replicate's
gate fluorescence; with evenly spaced gates assumed instead, every variant's K_D shifts by replicate
(+0.29, -0.12, -0.23 decades), which their fits do not show. The gate fluorescence is the mean log10
PE-A of the recorded cells inside each gate, under the authors' own gating, the definition of Phillips
2021's `fluor_mean_std_*.csv`.

Synonymous mode is on: each nucleotide variant (CDR1H + CDR3H nucleotides) is fitted on its own in each
replicate, and the synonymous variants of one amino-acid variant combine inside each replicate before the replicates
combine. This is the path the block takes on Adams through Amplicon Profiling.
"""

from __future__ import annotations

import csv
import math
import statistics
import urllib.parse
from functools import cache

import numpy as np
import polars as pl
from golden_support import cached_run, study_file
from paper_adams2016 import (
    CONCENTRATION_LABELS,
    CONCENTRATIONS_M,
    FSC_SSC_POLYGON,
    LOG10_KD_LOWER,
    LOG10_KD_UPPER,
    LOG10_PE_GATES,
    N_BINS,
    PRE_SORT_POLYGON,
    RAW_URL,
    SORT_COUNT_FILES,
    STUDY,
    TABLE_I,
    fcs_name,
)

from constants import FitParams
from replicate_combine import ReplicateInputs
from synonymous import SynonymousInputs

# The block's defaults for a new project (model/src/fitSettings.ts, kind/src/index.ts).
REPLICATES = ReplicateInputs()
PARAMS = FitParams(equation="log_fluorescence", plateau_floor=0.667, r2_gate_min=0.8, kd_se_gate_max=1.0)


def _file(name: str, folder: str) -> str:
    return str(study_file(STUDY, f"{folder}/{name}", RAW_URL.format(name=name)))


def _sort_fractions(rep: int) -> list[list[float]]:
    """Each gate's share of the cells sorted at each concentration, rows in CONCENTRATION_LABELS order."""
    text = open(_file(SORT_COUNT_FILES[rep], "extra")).read().replace("\r", "\n")
    rows = [ln.split("\t") for ln in text.split("\n") if ln.strip()][2:]
    by_label = {r[0].strip(): [float(x) for x in r[1 : 1 + N_BINS]] for r in rows}
    out = []
    for label in CONCENTRATION_LABELS:
        cells = by_label["0 M" if label == "0" else label.replace("10^", "")]
        out.append([c / sum(cells) for c in cells])
    return out


def _fcs_file(rep: int, concentration_index: int) -> str:
    name = fcs_name(rep, concentration_index)
    local = f"extra/fcs/rep{rep}_sort{concentration_index + 2}.fcs"
    return str(study_file(STUDY, local, RAW_URL.format(name=urllib.parse.quote(name))))


def _read_fcs(path: str) -> dict[str, np.ndarray]:
    """Channel name -> its value per recorded cell, from an FCS 3.0 file of 32-bit floats."""
    raw = open(path, "rb").read()
    text = raw[int(raw[10:18]) : int(raw[18:26]) + 1].decode("latin-1")
    fields = text[1:].split(text[0])
    meta = {k: v.strip() for k, v in zip(fields[::2], fields[1::2])}
    assert meta["$DATATYPE"] == "F" and meta["$MODE"] == "L", path
    n_par, n_cells = int(meta["$PAR"]), int(meta["$TOT"])
    assert all(meta[f"$P{i}B"] == "32" for i in range(1, n_par + 1)), path
    order = {"1,2,3,4": "<", "4,3,2,1": ">"}[meta["$BYTEORD"]]
    cells = np.frombuffer(raw, dtype=f"{order}f4", count=n_cells * n_par, offset=int(meta["$BEGINDATA"]))
    cells = cells.reshape(n_cells, n_par).astype(np.float64)
    return {meta[f"$P{i}N"]: cells[:, i - 1] for i in range(1, n_par + 1)}


def _inside(x: np.ndarray, y: np.ndarray, polygon: list[tuple[float, float]]) -> np.ndarray:
    """Which points lie inside a closed polygon: a ray to the right crosses its edges an odd number of times."""
    inside = np.zeros(x.shape, dtype=bool)
    for (x0, y0), (x1, y1) in zip(polygon, polygon[1:]):
        crosses = (y0 > y) != (y1 > y)
        at = np.where(crosses, x0 + (y - y0) * (x1 - x0) / (y1 - y0 if y1 != y0 else 1.0), -np.inf)
        inside ^= crosses & (x < at)
    return inside


@cache
def _gate_fluorescence(rep: int) -> list[list[float]]:
    """Each gate's mean log10 PE-A at each concentration, rows in CONCENTRATION_LABELS order.

    A gate that recorded no cell (gate 4 at 0 M in replicates 1 and 2) takes the replicate's mean
    position for that gate over the other concentrations: the block needs a position for every gate.
    """
    log10 = lambda v: np.log10(np.maximum(v, 1e-10))  # noqa: E731
    out: list[list[float | None]] = []
    for ci in range(len(CONCENTRATION_LABELS)):
        cells = _read_fcs(_fcs_file(rep, ci))
        keep = _inside(log10(cells["FSC-A"]), log10(cells["SSC-A"]), FSC_SSC_POLYGON)
        keep &= _inside(log10(cells["FITC-A"]), log10(cells["PE-A"]), PRE_SORT_POLYGON)
        pe = log10(cells["PE-A"][keep])
        gates = [pe[(pe > lo) & (pe < hi)] for lo, hi in LOG10_PE_GATES[rep]]
        out.append([float(g.mean()) if g.size else None for g in gates])
    for b in range(N_BINS):
        seen = [row[b] for row in out if row[b] is not None]
        for row in out:
            if row[b] is None:
                row[b] = statistics.fmean(seen)
    return out


@cache
def _tables() -> dict[int, pl.DataFrame]:
    return {
        rep: pl.read_csv(_file(f"replicate_{rep}.csv", "answer-key"), infer_schema_length=10000).with_columns(
            (pl.col("CDR1H") + "_" + pl.col("CDR3H")).alias("ntKey"),
            (pl.col("CDR1H_AA") + "_" + pl.col("CDR3H_AA")).alias("aaKey"),
        )
        for rep in (1, 2, 3)
    }


def reads() -> pl.DataFrame:
    parts = []
    for rep, table in _tables().items():
        fractions = _sort_fractions(rep)
        positions = _gate_fluorescence(rep)
        for ci, (label, conc) in enumerate(zip(CONCENTRATION_LABELS, CONCENTRATIONS_M)):
            for b in range(N_BINS):
                parts.append(
                    table.select(
                        pl.col("ntKey").alias("clonotypeKey"),
                        pl.lit(f"r{rep}_c{ci}_b{b + 1}").alias("sampleId"),
                        pl.lit(repr(conc) if conc else "0").alias("concentrationStr"),
                        pl.lit(conc, dtype=pl.Float64).alias("concentration"),
                        pl.lit(b + 1, dtype=pl.Int64).alias("bin"),
                        pl.col(f"fluorescein{label}bin{b}").fill_null(0).round(0).cast(pl.Int64).alias("reads"),
                        pl.lit(fractions[ci][b], dtype=pl.Float64).alias("sort_fraction"),
                        pl.lit(positions[ci][b], dtype=pl.Float64).alias("bin_fluorescence"),
                        pl.lit(str(rep)).alias("replicate"),
                    ).filter(pl.col("reads") > 0)
                )
    return pl.concat(parts)


def nt_to_aa() -> pl.DataFrame:
    return pl.concat([t.select("ntKey", "aaKey") for t in _tables().values()]).unique()


WILD_TYPE = "pJK36"
PARENT_ID = "wt"


def wild_type_links() -> tuple[pl.DataFrame, pl.DataFrame]:
    """Every amino-acid variant linked to the wild type, and its amino-acid distance from it.

    The study has one parent, so the link table is one parent id for every variant. The mutation
    count is the number of CDR1H and CDR3H positions that differ from the wild type's.
    """
    pairs = nt_to_aa()
    wild_type = pairs.filter(pl.col("ntKey") == clone_keys()[WILD_TYPE])["aaKey"].item()
    keys = sorted(pairs["aaKey"].unique().to_list())
    assert all(len(k) == len(wild_type) for k in keys)
    distance = [sum(a != b for a, b in zip(k, wild_type)) for k in keys]
    links = pl.DataFrame({"aaKey": keys, "parentId": [PARENT_ID] * len(keys)})
    mutations = pl.DataFrame({"aaKey": keys, "mutationCount": distance}, schema_overrides={"mutationCount": pl.Int64})
    return links, mutations


def run_block(params: FitParams = PARAMS, *, parent: bool = False) -> dict[str, pl.DataFrame]:
    """`parent` links every amino-acid variant to the wild type, so the run reports parent deltas."""
    from pipeline import run

    def go() -> dict:
        synonymous = SynonymousInputs(nt_to_aa=nt_to_aa())
        if parent:
            links, mutations = wild_type_links()
            synonymous = SynonymousInputs(nt_to_aa=nt_to_aa(), aa_parent_links=links, aa_mutations=mutations)
        return run(
            reads(),
            params=params,
            sort_fraction_col="sort_fraction",
            replicates=REPLICATES,
            synonymous=synonymous,
        )

    inputs = [_file(f"replicate_{r}.csv", "answer-key") for r in (1, 2, 3)]
    inputs += [_file(name, "extra") for name in SORT_COUNT_FILES.values()]
    inputs += [_fcs_file(r, ci) for r in (1, 2, 3) for ci in range(len(CONCENTRATION_LABELS))]
    inputs += [_file("titration_curves.csv", "extra")] if parent else []
    return cached_run("adams2016-parent" if parent else "adams2016", go, params, REPLICATES, inputs=inputs)


def _clip(x: float) -> float:
    return min(max(x, LOG10_KD_LOWER), LOG10_KD_UPPER)


def published() -> tuple[dict[str, float], dict[str, float]]:
    """The authors' per-replicate fits, combined into one value per variant.

    Per nucleotide variant, as their code does: each replicate's fit_KD clipped to the limits, then
    the mean log10 K_D over the three replicates.

    Per amino-acid variant the authors published no value: their Table I is not averaged over
    synonymous variants. This key is the test's own construction from their fits, in the order the
    block combines: the mean over a replicate's synonymous variants, then the mean over the
    replicates that have one. A difference from the block's amino-acid value therefore comes from
    the fits, not from the order of averaging.
    """
    acc_nt: dict[str, list[float]] = {}
    per_replicate_aa: list[dict[str, list[float]]] = []
    for table in _tables().values():
        acc_aa: dict[str, list[float]] = {}
        for key, akey, kd in table.select("ntKey", "aaKey", "fit_KD").iter_rows():
            if kd is None or not kd > 0 or not math.isfinite(kd):
                continue
            v = _clip(math.log10(kd))
            acc_nt.setdefault(key, []).append(v)
            acc_aa.setdefault(akey, []).append(v)
        per_replicate_aa.append(acc_aa)
    nt = {k: statistics.fmean(v) for k, v in acc_nt.items() if len(v) == 3}
    replicate_means: dict[str, list[float]] = {}
    for acc_aa in per_replicate_aa:
        for akey, values in acc_aa.items():
            replicate_means.setdefault(akey, []).append(statistics.fmean(values))
    aa = {k: statistics.fmean(v) for k, v in replicate_means.items()}
    return nt, aa


def clone_keys() -> dict[str, str]:
    """Table I clone name -> its nucleotide key, from data/titration_curves.csv."""
    path = _file("titration_curves.csv", "extra")
    out = {}
    for r in csv.DictReader(open(path)):
        r = {k.strip(): v.strip() for k, v in r.items()}
        out.setdefault(r["name"], f"{r['CDR1H']}_{r['CDR3H']}")
    return {c.name: out[c.name] for c in TABLE_I if c.name in out}


def in_range(x: float) -> bool:
    return LOG10_KD_LOWER < x < LOG10_KD_UPPER - 0.01
