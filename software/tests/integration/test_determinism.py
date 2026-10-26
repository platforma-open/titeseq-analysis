"""Identical inputs must produce byte-identical outputs across runs."""

from __future__ import annotations

import filecmp
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from helpers import cli_outputs, hill_reads

# Path to main.py — used to spawn the binary in a fresh subprocess so each run
# gets its own randomized ahash seed (the same condition that production runs
# meet, since each fit-curves invocation is a separate Python process).
_MAIN_PY = Path(__file__).resolve().parents[2] / "src" / "main.py"


_CORPUS_BIN_MODE = Path(__file__).resolve().parents[1] / "data" / "corpus" / "reads_bin_mode.parquet"


def _run_cli(reads_path: Path, out_dir: Path, extra_args: list[str]) -> dict[str, Path]:
    out_dir.mkdir()
    out_args, paths = cli_outputs(out_dir)
    result = subprocess.run(
        [sys.executable, str(_MAIN_PY), "--reads", str(reads_path), *out_args, *extra_args],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"fit-curves failed (rc={result.returncode}): {result.stderr}"
    return paths


def _default_model_library(tmp_path: Path) -> tuple[Path, list[str]]:
    # Needs a 0 M control. 410 clonotypes so the fit runs in parallel chunks, whose reassembly order must
    # not depend on which worker finishes first.
    concs = [0.0, 1e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 1e-6]
    rows = [
        row
        for i, kd in enumerate(np.logspace(-10, -6, 410))
        for row in hill_reads(f"CLONE_{i:03d}", float(kd), concs, [1, 2, 3, 4], cells=lambda b: 1000.0 + 100.0 * b)
    ]
    reads_path = tmp_path / "reads.parquet"
    pl.DataFrame(rows).write_parquet(reads_path)
    return reads_path, ["--allow-computed-baseline", "--cell-count-column", "cells"]


def _linear_hill_corpus(tmp_path: Path) -> tuple[Path, list[str]]:
    # The linear Hill path runs its own process pool over clonotypes.
    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps({"fit_model": "hill_free_n"}))
    return _CORPUS_BIN_MODE, ["--params", str(params_path)]


@pytest.mark.parametrize(
    "make_input", [_default_model_library, _linear_hill_corpus], ids=["default-model", "linear-hill"]
)
def test_outputs_are_byte_identical_across_runs(tmp_path, make_input):
    """Spawn fit-curves twice as separate processes and byte-compare every TSV.

    Each subprocess gets a fresh ahash seed, which is what surfaces hash-order
    non-determinism in `polars.group_by` — the actual condition under which
    CIDConflictError fires in production.
    """
    reads_path, extra_args = make_input(tmp_path)
    paths_a = _run_cli(reads_path, tmp_path / "run_a", extra_args)
    paths_b = _run_cli(reads_path, tmp_path / "run_b", extra_args)
    for name, path_a in paths_a.items():
        # shallow=False forces full content compare, not stat-based.
        assert filecmp.cmp(path_a, paths_b[name], shallow=False), (
            f"{name} differs between two runs of identical input — "
            f"non-determinism in pipeline will cause CIDConflictError downstream"
        )
