"""Adams 2016 as the block sees it: the authors' read counts per nucleotide variant, three replicates,
four gates with their sort fractions, run through `pipeline.run` with the block's default settings.

Synonymous mode is on: each nucleotide variant (CDR1H + CDR3H nucleotides) is fitted on its own in each
replicate, and the synonymous variants of one amino-acid variant combine inside each replicate before the replicates
combine. This is the path the block takes on Adams through Amplicon Profiling.
"""

from __future__ import annotations

import csv
import math
import statistics
from functools import cache

import polars as pl
from golden_support import cached_run, study_file
from paper_adams2016 import (
    CONCENTRATION_LABELS,
    CONCENTRATIONS_M,
    LOG10_KD_LOWER,
    LOG10_KD_UPPER,
    N_BINS,
    RAW_URL,
    SORT_COUNT_FILES,
    STUDY,
    TABLE_I,
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
                        pl.lit(str(rep)).alias("replicate"),
                    ).filter(pl.col("reads") > 0)
                )
    return pl.concat(parts)


def nt_to_aa() -> pl.DataFrame:
    return pl.concat([t.select("ntKey", "aaKey") for t in _tables().values()]).unique()


def run_block(params: FitParams = PARAMS) -> dict[str, pl.DataFrame]:
    from pipeline import run

    def go() -> dict:
        return run(
            reads(),
            params=params,
            sort_fraction_col="sort_fraction",
            replicates=REPLICATES,
            synonymous=SynonymousInputs(nt_to_aa=nt_to_aa()),
        )

    inputs = [_file(f"replicate_{r}.csv", "answer-key") for r in (1, 2, 3)]
    inputs += [_file(name, "extra") for name in SORT_COUNT_FILES.values()]
    return cached_run("adams2016", go, params, REPLICATES, inputs=inputs)


def _clip(x: float) -> float:
    return min(max(x, LOG10_KD_LOWER), LOG10_KD_UPPER)


def published() -> tuple[dict[str, float], dict[str, float]]:
    """The authors' values, combined as their code does: each replicate's fit_KD clipped to the limits,
    then the mean log10 K_D. Per nucleotide variant (over replicates) and per amino-acid variant (over
    replicates and synonymous variants)."""
    nt, aa = {}, {}
    acc_nt: dict[str, list[float]] = {}
    acc_aa: dict[str, list[float]] = {}
    for table in _tables().values():
        for key, akey, kd in table.select("ntKey", "aaKey", "fit_KD").iter_rows():
            if kd is None or not kd > 0 or not math.isfinite(kd):
                continue
            v = _clip(math.log10(kd))
            acc_nt.setdefault(key, []).append(v)
            acc_aa.setdefault(akey, []).append(v)
    nt = {k: statistics.fmean(v) for k, v in acc_nt.items() if len(v) == 3}
    aa = {k: statistics.fmean(v) for k, v in acc_aa.items()}
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
