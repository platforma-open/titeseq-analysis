"""Identical inputs must produce byte-identical outputs across runs."""

from __future__ import annotations

import filecmp
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
from helpers import hill_reads

# Path to main.py — used to spawn the binary in a fresh subprocess so each run
# gets its own randomized ahash seed (the same condition that production runs
# meet, since each fit-curves invocation is a separate Python process).
_MAIN_PY = Path(__file__).resolve().parents[2] / "src" / "main.py"


def _run_cli(reads_path: Path, out_dir: Path) -> dict[str, Path]:
    pc = out_dir / "per_clonotype.tsv"
    mb = out_dir / "mean_bin.tsv"
    fmb = out_dir / "fitted_mean_bin.tsv"
    cv = out_dir / "concentration_value.tsv"
    result = subprocess.run(
        [
            sys.executable,
            str(_MAIN_PY),
            "--reads",
            str(reads_path),
            "--out-per-clonotype",
            str(pc),
            "--out-mean-bin",
            str(mb),
            "--out-fitted-mean-bin",
            str(fmb),
            "--out-concentration-value",
            str(cv),
            "--cell-count-column",
            "cells",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"fit-curves failed (rc={result.returncode}): {result.stderr}"
    return {"per_clonotype": pc, "mean_bin": mb, "fitted_mean_bin": fmb, "concentration_value": cv}


def test_outputs_are_byte_identical_across_runs(tmp_path):
    """Spawn fit-curves twice as separate processes and byte-compare every TSV.

    Each subprocess gets a fresh ahash seed, which is what surfaces hash-order
    non-determinism in `polars.group_by` — the actual condition under which
    CIDConflictError fires in production.
    """
    # Default model (needs a 0 M control). 410 clonotypes so the fit runs in parallel chunks, whose
    # reassembly order must not depend on which worker finishes first.
    concs = [0.0, 1e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 1e-6]
    bins = [1, 2, 3, 4]
    rows: list[dict] = []
    for i, kd in enumerate(np.logspace(-10, -6, 410)):
        rows.extend(hill_reads(f"CLONE_{i:03d}", float(kd), concs, bins, cells=lambda b: 1000.0 + 100.0 * b))

    reads_path = tmp_path / "reads.parquet"
    pl.DataFrame(rows).write_parquet(reads_path)

    out_a = tmp_path / "run_a"
    out_b = tmp_path / "run_b"
    out_a.mkdir()
    out_b.mkdir()

    paths_a = _run_cli(reads_path, out_a)
    paths_b = _run_cli(reads_path, out_b)

    for name, path_a in paths_a.items():
        path_b = paths_b[name]
        assert path_a.exists(), f"first run missing {name}"
        assert path_b.exists(), f"second run missing {name}"
        # shallow=False forces full content compare, not stat-based.
        assert filecmp.cmp(path_a, path_b, shallow=False), (
            f"{name} differs between two runs of identical input — "
            f"non-determinism in pipeline will cause CIDConflictError downstream"
        )
