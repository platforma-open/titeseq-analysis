"""Adams 2016 golden data: download the authors' processed tables at a pinned commit and build the
block's reads frame from them.

The files are cached under $TITESEQ_GOLDEN_CACHE (default ~/.cache/titeseq-analysis/adams2016). If they
cannot be fetched, every golden test is skipped with the reason; nothing from the download is committed.
"""

from __future__ import annotations

import io
import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pytest
from paper_adams2016 import CONCENTRATION_LABELS, CONCENTRATIONS_M, N_BINS, PINNED_COMMIT, REPLICATES, REPO

RAW_URL = f"https://raw.githubusercontent.com/{REPO}/{PINNED_COMMIT}/data/{{name}}"


def _cache_dir() -> Path:
    root = os.environ.get("TITESEQ_GOLDEN_CACHE", str(Path.home() / ".cache" / "titeseq-analysis" / "adams2016"))
    path = Path(root) / PINNED_COMMIT[:12]
    path.mkdir(parents=True, exist_ok=True)
    return path


def _fetch(name: str) -> Path:
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
    lines = [ln for ln in _fetch(name).read_text().replace("\r", "\n").split("\n") if ln.strip()]
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


def _build(rep) -> GoldenReplicate:
    table = pl.read_csv(_fetch(f"replicate_{rep.number}.csv"), infer_schema_length=10000)
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
    return GoldenReplicate(rep.number, pl.concat(parts), published, rep.bin_width_decades, rep.min_plateau_fraction)


@pytest.fixture(scope="session")
def adams_replicates() -> list[GoldenReplicate]:
    return [_build(rep) for rep in REPLICATES]


@pytest.fixture(scope="session")
def adams_clone_keys() -> dict[str, str]:
    """Table I clone name -> clonotypeKey (CDR1H_CDR3H nucleotides), from the authors' flow data."""
    # The file uses classic-Mac "\r" line endings.
    text = _fetch("titration_curves.csv").read_text().replace("\r\n", "\n").replace("\r", "\n")
    curves = pl.read_csv(io.StringIO(text), infer_schema=False)
    curves = curves.rename({c: c.strip() for c in curves.columns})
    return {
        str(row["name"]).strip(): f"{str(row['CDR1H']).strip()}_{str(row['CDR3H']).strip()}"
        for row in curves.unique(subset="name", keep="first").iter_rows(named=True)
    }
