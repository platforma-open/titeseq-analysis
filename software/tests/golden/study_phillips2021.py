"""Phillips 2021 (CR6261) as the block sees it: the authors' read counts per 11-bit genotype, one campaign
(H1 or H9) at a time, its three replicate libraries as sample replicates, 12 concentrations (0 M included)
x 4 gates, with each gate's sort fraction and measured log10 fluorescence bound. One sequence per
genotype, so no synonymous mode.

Inputs:
  reads             extra/author-counts/df_merged_6261<lib>.csv, `correct_counts_6261_<lib>_<code>_<bin>`
                    (fractional after the authors' correction; rounded to whole reads)
  sort fraction,    upload/metadata.csv, the per-sample metadata as the block receives it. The loader checks
  gate fluorescence it against the authors' own tables (extra/cell_counts/Binding: cells per gate;
                    extra/gate-fluorescence: mean log10 fluorescence per gate) and fails on any mismatch.
  concentration     the metadata's labels (10^-7 ... 10^-12 M), or the authors' fitting values 0.05 decades
                    lower (`shifted=True`). Scaling every concentration by one factor scales every fitted
                    K_D by that factor, so the shift moves each log10 K_D by exactly 0.05.
"""

from __future__ import annotations

import csv
import math
import statistics
from functools import cache

import polars as pl
from golden_support import cached_run, study_file
from paper_phillips2021 import (
    AUTHOR_CONC_SHIFT,
    LIBRARIES,
    N_BINS,
    NEG_LOG10_CONC,
    NEG_LOG10_KD_LOWER,
    NEG_LOG10_KD_UPPER,
    STUDY,
)

from constants import FitParams
from replicate_combine import ReplicateInputs

REPLICATES = ReplicateInputs()
# The block's defaults for a new project (model/src/fitSettings.ts, kind/src/index.ts).
PARAMS = FitParams(equation="log_fluorescence", plateau_floor=0.667, r2_gate_min=0.8, kd_se_gate_max=1.0)


def _file(path: str) -> str:
    return str(study_file(STUDY, path))


def _code_of(conc_m: str) -> str:
    """The authors' concentration code for a metadata molar value: '3.162e-08' -> '75'."""
    value = float(conc_m)
    if value == 0:
        return "0"
    neg = round(-math.log10(value), 2)
    return next(code for code, x in NEG_LOG10_CONC.items() if x is not None and abs(x - neg) < 1e-3)


@cache
def _metadata() -> dict[tuple[str, str, int], dict]:
    """(library, concentration code, gate) -> the sample's metadata row."""
    out = {}
    for r in csv.DictReader(open(_file("upload/metadata.csv"))):
        out[(r["Replicate"], _code_of(r["Antigen concentration (M)"]), int(r["FACS bin"]))] = r
    return out


def _authors_gates(lib: str) -> dict[str, tuple[list[float], list[float]]]:
    """code -> (sort fractions, gate log10 fluorescence) from the authors' own tables."""
    cells = {}
    for row in list(csv.reader(open(_file(f"extra/cell_counts/Binding/6261{lib}_cellct.csv"))))[1:]:
        counts = [float(x) for x in row[1 : 1 + N_BINS]]
        cells[row[0].strip("_")] = [c / sum(counts) for c in counts]
    fluor = {}
    for row in csv.reader(open(_file(f"extra/gate-fluorescence/fluor_mean_std_6261{lib}.csv"))):
        fluor[row[0].strip("_")] = [float(x) for x in row[1 : 1 + N_BINS]]
    return {code: (cells[code], fluor[code]) for code in NEG_LOG10_CONC}


def _check_metadata(lib: str) -> None:
    meta = _metadata()
    for code, (fractions, fluor) in _authors_gates(lib).items():
        for b in range(N_BINS):
            m = meta[(lib, code, b + 1)]
            assert math.isclose(float(m["Sort fraction"]), fractions[b], abs_tol=1e-9), (lib, code, b)
            assert math.isclose(float(m["Gate fluorescence (log10)"]), fluor[b], abs_tol=1e-9), (lib, code, b)


def concentration(code: str, shifted: bool) -> tuple[str, float]:
    """(concentrationStr, molar value) for one concentration code."""
    neg = NEG_LOG10_CONC[code]
    if neg is None:
        return "0", 0.0
    if not shifted:
        m = next(r for r in _metadata().values() if _code_of(r["Antigen concentration (M)"]) == code)
        return m["Antigen concentration (M)"], float(m["Antigen concentration (M)"])
    value = 10 ** -(neg + AUTHOR_CONC_SHIFT)
    return repr(value), value


@cache
def _counts(lib: str) -> pl.DataFrame:
    return pl.read_csv(
        _file(f"extra/author-counts/df_merged_6261{lib}.csv"),
        schema_overrides={"variant": pl.Utf8},
        infer_schema_length=10000,
    )


@cache
def reads(campaign: str, shifted: bool) -> pl.DataFrame:
    meta = _metadata()
    parts = []
    for lib in LIBRARIES[campaign]:
        _check_metadata(lib)
        table = _counts(lib)
        assert table["variant"].str.len_chars().unique().to_list() == [11], lib
        for code in NEG_LOG10_CONC:
            conc_str, conc = concentration(code, shifted)
            for b in range(1, N_BINS + 1):
                m = meta[(lib, code, b)]
                parts.append(
                    table.select(
                        pl.col("variant").alias("clonotypeKey"),
                        pl.lit(m["Sample"]).alias("sampleId"),
                        pl.lit(conc_str).alias("concentrationStr"),
                        pl.lit(conc, dtype=pl.Float64).alias("concentration"),
                        pl.lit(b, dtype=pl.Int64).alias("bin"),
                        pl.col(f"correct_counts_6261_{lib}_{code}_{b}")
                        .fill_null(0)
                        .round(0)
                        .cast(pl.Int64)
                        .alias("reads"),
                        pl.lit(float(m["Sort fraction"]), dtype=pl.Float64).alias("sort_fraction"),
                        pl.lit(float(m["Gate fluorescence (log10)"]), dtype=pl.Float64).alias("bin_fluorescence"),
                        pl.lit(lib).alias("replicate"),
                    ).filter(pl.col("reads") > 0)
                )
    return pl.concat(parts)


def run_block(campaign: str, params: FitParams = PARAMS, shifted: bool = False) -> dict[str, pl.DataFrame]:
    from pipeline import run

    def go() -> dict:
        return run(reads(campaign, shifted), params=params, sort_fraction_col="sort_fraction", replicates=REPLICATES)

    inputs = [_file("upload/metadata.csv")]
    for lib in LIBRARIES[campaign]:
        inputs += [
            _file(f"extra/author-counts/df_merged_6261{lib}.csv"),
            _file(f"extra/cell_counts/Binding/6261{lib}_cellct.csv"),
            _file(f"extra/gate-fluorescence/fluor_mean_std_6261{lib}.csv"),
        ]
    return cached_run(f"phillips2021-{campaign}", go, params, REPLICATES, shifted, inputs=inputs)


@cache
def published(campaign: str) -> dict[str, float]:
    """genotype -> the authors' unfiltered mean -log10 K_D over the replicates
    (answer-key/20210323_<campaign>_avg_unf.csv: the plain mean of the replicates that have a value)."""
    out = {}
    for r in csv.DictReader(open(_file(f"answer-key/20210323_{campaign}_avg_unf.csv"))):
        v = r[f"{campaign}_mean"]
        if v not in ("", "nan"):
            out[r["variant"]] = float(v)
    return out


@cache
def isogenic(campaign: str) -> dict[str, float]:
    """Isogenic clone name -> its flow-cytometry log10 K_D: the mean of its two replicate titrations."""
    acc: dict[str, list[float]] = {}
    for r in csv.DictReader(open(_file("answer-key/isogenic/Isogenic_Kd/Mutants/fitted_Kds.csv"))):
        if r["antigen"].lower() == campaign:
            acc.setdefault(r["shortname"], []).append(float(r["KD"]))
    return {k: statistics.fmean(v) for k, v in acc.items()}


def clip(neg_log10_kd: float) -> float:
    return min(max(neg_log10_kd, NEG_LOG10_KD_LOWER), NEG_LOG10_KD_UPPER)


def in_range(neg_log10_kd: float) -> bool:
    return NEG_LOG10_KD_LOWER + 0.01 < neg_log10_kd < NEG_LOG10_KD_UPPER - 0.01
