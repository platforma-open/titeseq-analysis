"""Shared machinery for the published-study golden tests.

Study data carries no licence, so nothing from it is committed. Each file is looked up, in order:

  1. $TITESEQ_STUDIES_DIR/<study>/<path>
  2. a `.data/titeseq-studies/<study>/<path>` in any folder above this one (the mictx workspace)
  3. $TITESEQ_GOLDEN_CACHE (default ~/.cache/titeseq-analysis/golden)/<study>/<path>
  4. a download from the authors' public repository at a pinned commit, into that cache

and the test is skipped, with the reason, when none of them yields the file.

What is committed is in `pins/<case>.json`: agreement with the published values, the per-clone values
the papers measured one by one, and a digest of the whole output. The whole output itself is kept in the
cache, keyed by its digest. When the digest moves, the test fails with a diff against that snapshot.
`TITESEQ_GOLDEN_UPDATE=1` rewrites the pins from the current code instead of comparing.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import statistics
import urllib.error
import urllib.request
from collections import Counter
from collections.abc import Callable
from pathlib import Path

import numpy as np
import polars as pl
import pytest
import scipy

HERE = Path(__file__).resolve().parent
PINS_DIR = HERE / "pins"
SRC_DIR = HERE.parents[1] / "src"
UPDATE = os.environ.get("TITESEQ_GOLDEN_UPDATE") == "1"


def cache_root() -> Path:
    root = Path(os.environ.get("TITESEQ_GOLDEN_CACHE", str(Path.home() / ".cache" / "titeseq-analysis" / "golden")))
    root.mkdir(parents=True, exist_ok=True)
    return root


def _workspace_study_dir(study: str) -> Path | None:
    for parent in HERE.parents:
        candidate = parent / ".data" / "titeseq-studies" / study
        if candidate.is_dir():
            return candidate
    return None


def study_file(study: str, path: str, url: str | None = None) -> Path:
    """The local copy of one study file; downloaded into the cache when `url` is given and no copy exists."""
    roots = []
    if os.environ.get("TITESEQ_STUDIES_DIR"):
        roots.append(Path(os.environ["TITESEQ_STUDIES_DIR"]) / study)
    workspace = _workspace_study_dir(study)
    if workspace is not None:
        roots.append(workspace)
    cached = cache_root() / study
    roots.append(cached)
    for root in roots:
        if (root / path).is_file():
            return root / path
    if url is None:
        pytest.skip(f"{study} data not found: {path} (set TITESEQ_STUDIES_DIR)")
    target = cached / path
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            data = resp.read()
    except (urllib.error.URLError, TimeoutError) as exc:
        pytest.skip(f"{study} data unavailable ({path}): {exc}")
    tmp = target.with_suffix(target.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(target)
    return target


def _code_key(*parts: object) -> str:
    """Everything that can change a run: the software's source, the library versions, the inputs."""
    h = hashlib.sha256()
    for f in sorted(SRC_DIR.glob("*.py")):
        h.update(f.name.encode() + f.read_bytes())
    for lib in (np, scipy, pl):
        h.update(f"{lib.__name__}={lib.__version__}".encode())
    h.update(repr(parts).encode())
    return h.hexdigest()[:24]


# Bumped when the cached layout changes, so older cache folders are not read.
CACHE_LAYOUT = 2


def cached_run(case: str, run: Callable[[], dict], *key_parts: object) -> dict:
    """`run()`'s frame and dict outputs, cached on disk until the source, a library or `key_parts` change.
    Frames are kept as parquet, dicts as JSON; every other output is dropped."""
    folder = cache_root() / "runs" / f"{case}-{_code_key(case, CACHE_LAYOUT, *key_parts)}"
    if folder.is_dir() and (folder / "_done").exists():
        outputs: dict = {p.stem: pl.read_parquet(p) for p in folder.glob("*.parquet")}
        outputs.update({p.stem: json.loads(p.read_text()) for p in folder.glob("*.json")})
        return outputs
    outputs = {k: v for k, v in run().items() if isinstance(v, (pl.DataFrame, dict))}
    folder.mkdir(parents=True, exist_ok=True)
    for name, value in outputs.items():
        if isinstance(value, pl.DataFrame):
            value.write_parquet(folder / f"{name}.parquet")
        else:
            (folder / f"{name}.json").write_text(json.dumps(value, sort_keys=True))
    (folder / "_done").write_text("")
    return outputs


# ---- The snapshot: one row per reported variant, the columns a user reads.

SNAPSHOT_SCHEMA = {"key": pl.Utf8, "log10Kd": pl.Float64, "log10KdSe": pl.Float64, "affinityClass": pl.Utf8}


def snapshot(frame: pl.DataFrame, key: str) -> pl.DataFrame:
    """[key, log10Kd, log10KdSe, affinityClass], sorted by key, values rounded to 1e-6 decades."""
    return (
        frame.select(
            pl.col(key).cast(pl.Utf8).alias("key"),
            pl.col("kd").log10().round(6).alias("log10Kd"),
            pl.col("kdLog10Se").round(6).alias("log10KdSe"),
            pl.col("affinityClass"),
        )
        .cast(SNAPSHOT_SCHEMA)
        .sort("key")
    )


def digest(snap: pl.DataFrame) -> str:
    h = hashlib.sha256()
    for row in snap.iter_rows():
        h.update(repr(row).encode())
    return h.hexdigest()[:16]


def _snapshot_path(case: str, dig: str) -> Path:
    return cache_root() / "snapshots" / f"{case}-{dig}.parquet"


def diff_report(before: pl.DataFrame, after: pl.DataFrame, limit: int = 12) -> str:
    """What moved between two snapshots: class moves, Kd,app and SE changes, the largest of them."""
    joined = before.join(after, on="key", how="full", suffix="_new", coalesce=True)
    lines = ["pinned -> now", f"rows: {before.height} -> {after.height}"]
    moves = Counter((a, b) for a, b in joined.select("affinityClass", "affinityClass_new").iter_rows() if a != b)
    lines.append("class moves: " + (", ".join(f"{a} -> {b}: {n}" for (a, b), n in moves.most_common()) or "none"))
    delta = joined.with_columns(
        (pl.col("log10Kd_new") - pl.col("log10Kd")).alias("dKd"),
        (pl.col("log10KdSe_new") - pl.col("log10KdSe")).alias("dSe"),
    )
    kd_moved = delta.filter(pl.col("dKd").abs() > 1e-6)
    se_moved = delta.filter(pl.col("dSe").abs() > 1e-6)
    lines.append(f"Kd,app changed on {kd_moved.height} variants; SE changed on {se_moved.height}")
    for col in (c for c in before.columns if c in after.columns and c not in SNAPSHOT_SCHEMA):
        changed = joined.filter(~pl.col(col).eq_missing(pl.col(f"{col}_new"))).height
        if changed:
            lines.append(f"{col} changed on {changed} rows")
    if kd_moved.height:
        lines.append(f"  |dlog10Kd| median {kd_moved['dKd'].abs().median():.3f}, max {kd_moved['dKd'].abs().max():.3f}")
    top = delta.with_columns(pl.col("dKd").abs().fill_null(0).alias("_a")).sort("_a", descending=True).head(limit)
    for r in top.filter(pl.col("_a") > 1e-6).iter_rows(named=True):
        lines.append(
            f"  {r['key']}: {r['log10Kd']} -> {r['log10Kd_new']} ({r['affinityClass']} -> {r['affinityClass_new']})"
        )
    return "\n".join(lines)


# ---- Agreement with the published values.


def agreement(pairs: list[tuple[float, float]]) -> dict:
    """n, Pearson r and the median absolute difference of (ours, published) log10 K_D pairs."""
    if len(pairs) < 3:
        return {"n": len(pairs)}
    ours, pub = zip(*pairs)
    d = sorted(abs(a - b) for a, b in pairs)
    return {
        "n": len(pairs),
        "r": round(statistics.correlation(ours, pub), 4),
        "medianAbsDiff": round(d[len(d) // 2], 4),
        "within0.5": round(sum(v < 0.5 for v in d) / len(d), 4),
    }


# ---- Pins.


def load_pins(case: str) -> dict | None:
    path = PINS_DIR / f"{case}.json"
    return json.loads(path.read_text()) if path.exists() else None


def write_pins(case: str, pins: dict) -> None:
    PINS_DIR.mkdir(exist_ok=True)
    (PINS_DIR / f"{case}.json").write_text(json.dumps(pins, indent=2, sort_keys=True) + "\n")


def check_digest(case: str, snap: pl.DataFrame, pins: dict) -> None:
    """Fail with a diff when the whole output moved since the pins were written."""
    dig = digest(snap)
    path = _snapshot_path(case, dig)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        snap.write_parquet(path)
    if dig == pins["digest"]:
        return
    old = _snapshot_path(case, pins["digest"])
    report = (
        diff_report(pl.read_parquet(old), snap)
        if old.exists()
        else "no local snapshot of the pinned output to diff against: run the pinned commit once to create it"
    )
    pytest.fail(
        f"{case}: the output moved (digest {pins['digest']} -> {dig}).\n{report}\n"
        f"If the change is intended, rerun with TITESEQ_GOLDEN_UPDATE=1 and commit pins/{case}.json.",
        pytrace=False,
    )


def update_or_load(case: str, current: dict, snap: pl.DataFrame) -> dict:
    """Under TITESEQ_GOLDEN_UPDATE=1 write `current` as the pins and return it; else return the committed pins."""
    if UPDATE:
        path = _snapshot_path(case, current["digest"])
        path.parent.mkdir(parents=True, exist_ok=True)
        snap.write_parquet(path)
        write_pins(case, current)
        return current
    pins = load_pins(case)
    if pins is None:
        pytest.fail(f"no pins/{case}.json: create it with TITESEQ_GOLDEN_UPDATE=1", pytrace=False)
    return pins


def close(a: float | None, b: float | None, tol: float = 1e-6) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, abs_tol=tol)
