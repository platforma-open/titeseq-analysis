"""Adams 2016 golden data: download the authors' processed tables at a pinned commit, build the block's
reads frame from them, and run cached fits.

The files are cached under $TITESEQ_GOLDEN_CACHE (default ~/.cache/titeseq-analysis/adams2016). If they
cannot be fetched, every golden test is skipped with the reason; nothing from the download is committed.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from paper_adams2016 import (
    CONCENTRATION_LABELS,
    CONCENTRATIONS_M,
    N_BINS,
    PINNED_COMMIT,
    REPO,
    WT_CDR1H_AA,
    WT_CDR3H_AA,
)

RAW_URL = f"https://raw.githubusercontent.com/{REPO}/{PINNED_COMMIT}/data/{{name}}"


def _cache_dir() -> Path:
    root = os.environ.get("TITESEQ_GOLDEN_CACHE", str(Path.home() / ".cache" / "titeseq-analysis" / "adams2016"))
    path = Path(root) / PINNED_COMMIT[:12]
    path.mkdir(parents=True, exist_ok=True)
    return path


def fetch(name: str) -> Path:
    target = _cache_dir() / name
    if not target.exists():
        try:
            with urllib.request.urlopen(RAW_URL.format(name=name), timeout=60) as resp:
                data = resp.read()
        except (urllib.error.URLError, TimeoutError) as exc:
            pytest.skip(f"Adams 2016 data unavailable ({name}): {exc}")
        tmp = target.with_suffix(".part")
        tmp.write_bytes(data)
        tmp.replace(target)
    return target


def _sort_counts(name: str) -> list[list[float]]:
    """Cells sorted per (concentration, bin). Rows follow CONCENTRATION_LABELS; the expression row is dropped."""
    lines = [ln for ln in fetch(name).read_text().replace("\r", "\n").split("\n") if ln.strip()]
    rows = [ln.split("\t") for ln in lines[2:]]  # skip "#replicate N" and the header
    by_label = {r[0].strip(): [float(x) for x in r[1 : 1 + N_BINS]] for r in rows}
    out = [by_label["0 M"]]
    out += [by_label[label.replace("10^", "")] for label in CONCENTRATION_LABELS[1:]]
    return out


@dataclass(frozen=True)
class GoldenReplicate:
    number: int
    reads: pl.DataFrame  # block input: clonotypeKey, sampleId, concentrationStr, concentration, bin, reads, cell_count
    published: pl.DataFrame  # clonotypeKey, cdr1h_aa, cdr3h_aa, published_log10_kd (clipped to the limits)
    bin_width_decades: float
    min_plateau_fraction: float
    baseline_fluorescence: float


def build_replicate(rep) -> GoldenReplicate:
    table = pl.read_csv(fetch(f"replicate_{rep.number}.csv"), infer_schema_length=10000)
    table = table.with_columns((pl.col("CDR1H") + "_" + pl.col("CDR3H")).alias("clonotypeKey"))
    cells = _sort_counts(rep.sort_counts_file)
    parts = []
    for ci, (label, conc) in enumerate(zip(CONCENTRATION_LABELS, CONCENTRATIONS_M)):
        for b in range(N_BINS):
            parts.append(
                table.select(
                    "clonotypeKey",
                    pl.lit(f"c{ci}_b{b + 1}").alias("sampleId"),
                    pl.lit(repr(conc) if conc else "0").alias("concentrationStr"),
                    pl.lit(conc, dtype=pl.Float64).alias("concentration"),
                    pl.lit(b + 1, dtype=pl.Int64).alias("bin"),
                    pl.col(f"fluorescein{label}bin{b}").fill_null(0).round(0).cast(pl.Int64).alias("reads"),
                    pl.lit(cells[ci][b], dtype=pl.Float64).alias("cell_count"),
                ).filter(pl.col("reads") > 0)
            )
    published = table.select(
        "clonotypeKey",
        pl.col("CDR1H_AA").alias("cdr1h_aa"),
        pl.col("CDR3H_AA").alias("cdr3h_aa"),
        pl.col("fit_KD").log10().clip(-9.5, -5.0).alias("published_log10_kd"),
    )
    return GoldenReplicate(
        rep.number,
        pl.concat(parts),
        published,
        rep.bin_width_decades,
        rep.min_plateau_fraction,
        rep.baseline_fluorescence,
    )


_SRC = Path(__file__).resolve().parents[2] / "src"
REST = "__rest_of_library__"


def _fit_key(rep_number: int, params, keys: tuple[str, ...] | None) -> str:
    """Everything that can change a fit: source code, library versions, settings, data, subset."""
    h = hashlib.sha256()
    for f in sorted(_SRC.glob("*.py")):
        h.update(f.name.encode() + f.read_bytes())
    import scipy

    for lib in (np, scipy, pl):
        h.update(f"{lib.__name__}={lib.__version__}".encode())
    h.update(repr((PINNED_COMMIT, rep_number, params, keys)).encode())
    return h.hexdigest()[:24]


def fit_golden(rep: GoldenReplicate, params, keys: tuple[str, ...] | None = None) -> pl.DataFrame:
    """Run pipeline.run on a replicate and join the published K_D; cached on disk by `_fit_key`.

    `keys` restricts the fit to those clonotypes: every other clonotype is pooled into one, so each
    sample's total reads (and the bin exposures built from them) are unchanged.
    """
    from pipeline import run

    path = _cache_dir() / "fits" / f"rep{rep.number}-{_fit_key(rep.number, params, keys)}.parquet"
    if path.exists():
        return pl.read_parquet(path)
    reads = rep.reads
    if keys is not None:
        inside = pl.col("clonotypeKey").is_in(list(keys))
        rest = (
            reads.filter(~inside)
            .group_by("sampleId", "concentrationStr", "concentration", "bin", "cell_count")
            .agg(pl.col("reads").sum())
            .with_columns(pl.lit(REST).alias("clonotypeKey"))
            .select(reads.columns)
        )
        reads = pl.concat([reads.filter(inside), rest])
    with contextlib.redirect_stdout(io.StringIO()):
        out = run(reads, params=params, cell_count_col="cell_count")
    pc = out["per_clonotype"].with_columns(pl.col("kd").log10().alias("log10_kd"))
    fits = rep.published.join(pc, on="clonotypeKey", how="inner" if keys is not None else "left")
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.write_parquet(path)
    return fits


def fast_subset(rep: GoldenReplicate, clone_keys: dict[str, str], per_stratum: int = 30) -> tuple[str, ...]:
    """Table I clones, WT synonymous variants, and a seeded sample from each published-K_D stratum."""
    pub = rep.published.filter(pl.col("published_log10_kd").is_not_nan())
    keys = set(clone_keys.values())
    keys |= set(pub.filter((pl.col("cdr1h_aa") == WT_CDR1H_AA) & (pl.col("cdr3h_aa") == WT_CDR3H_AA))["clonotypeKey"])
    edges = [-10.0, -9.4, -8.5, -7.5, -6.5, -5.5, -5.01, 0.0]
    rng = np.random.default_rng(rep.number)
    for lo, hi in zip(edges[:-1], edges[1:]):
        stratum = sorted(pub.filter(pl.col("published_log10_kd").is_between(lo, hi, closed="right"))["clonotypeKey"])
        keys |= set(rng.choice(stratum, min(per_stratum, len(stratum)), replace=False).tolist())
    return tuple(sorted(keys & set(pub["clonotypeKey"])))
