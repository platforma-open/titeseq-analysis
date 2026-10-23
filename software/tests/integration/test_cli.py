"""CLI smoke tests for main.main(argv) — exercise argparse, file I/O, output writers."""

from __future__ import annotations

import json

import numpy as np
import polars as pl
import pytest
from helpers import cli_outputs, hill_reads

from main import main


def _cli_reads(clonotype: str, true_kd: float, concs: list[float], bins: list[int]) -> list[dict]:
    return hill_reads(clonotype, true_kd, concs, bins, cells=lambda b: 1000.0)


# Default model: needs a 0 M control and the cells sorted per sample.
_DEFAULT_MODEL_CONCS = [0.0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
_CELLS_ARGS = ["--cell-count-column", "cells"]


# Both writer branches: every output file is written in the format its suffix names.
@pytest.mark.parametrize("suffix", [".parquet", ".tsv"])
def test_cli_writes_every_output(tmp_path, suffix):
    reads_path = tmp_path / "reads.parquet"
    pl.DataFrame(_cli_reads("G1", 1e-8, _DEFAULT_MODEL_CONCS, [1, 2, 3, 4])).write_parquet(reads_path)
    out_args, paths = cli_outputs(tmp_path, suffix)
    assert main(["--reads", str(reads_path), *out_args, *_CELLS_ARGS]) == 0
    assert all(path.exists() for path in paths.values())
    pc_path = paths["per_clonotype"]
    pc = pl.read_parquet(pc_path) if suffix == ".parquet" else pl.read_csv(pc_path, separator="\t")
    assert pc["clonotypeKey"].to_list() == ["G1"]


# --params JSON — ensures FitParams override path is wired through argparse
def test_cli_with_params_json(tmp_path):
    rows = _cli_reads("G1", 1e-8, [1e-10, 1e-9, 1e-8, 1e-7, 1e-6], [1, 2, 3, 4])
    reads_path = tmp_path / "reads.parquet"
    pl.DataFrame(rows).write_parquet(reads_path)

    params_path = tmp_path / "params.json"
    params_path.write_text(
        json.dumps(
            {
                "min_reads_per_concentration": 3,
                "min_concentration_points": 5,
                "r2_threshold_good": 0.8,
                "r2_threshold_failed": 0.5,
                "n_min": 0.5,
                "n_max": 2.0,
                "hook_effect_threshold_bin": 0.2,
                "hook_effect_threshold_no_bin": 0.02,
                "hook_effect_min_reads": 20,
                "fit_model": "hill_free_n",
            }
        )
    )

    out_args, _ = cli_outputs(tmp_path)
    assert main(["--reads", str(reads_path), *out_args, "--params", str(params_path)]) == 0


# Hook-flagged clonotype → non_monotonic_signal branch in run()
def test_cli_hook_effect_triggers_failure(tmp_path):
    # R9b requires BOTH clauses: (top2-top1 > theta) AND (top3-top1 > theta/2).
    # Signals must stay elevated at top-3 AND top-2 before dropping at top-1.
    # Sub-µM grid spanning 0.1 nM → 1 µM.
    concs = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
    bins = [1, 2, 3, 4]
    rows = []
    # Top-3 (c=1e-8) and top-2 (c=1e-7) both populate bins {3,4} → signal 3.5;
    # top-1 (c=1e-6) populates bins {2,3} → signal 2.5 → drop of 1.0.
    targets = [1.5, 2.0, 3.5, 3.6, 2.5]
    for i, c in enumerate(concs):
        target = targets[i]
        weights = np.exp(-0.5 * ((np.array(bins, dtype=float) - target) / 0.25) ** 2)
        weights /= weights.sum()
        counts = np.round(weights * 500).astype(int)
        for j, b in enumerate(bins):
            rows.append(
                {
                    "clonotypeKey": "H1",
                    "sampleId": f"s_c{i}_b{b}",
                    "concentrationStr": str(c),
                    "concentration": float(c),
                    "bin": int(b),
                    "reads": int(counts[j]),
                }
            )
    reads_path = tmp_path / "reads.parquet"
    pl.DataFrame(rows).write_parquet(reads_path)

    # The hook gate belongs to the opt-in linear Hill model; the default model does not fail on it.
    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps({"fit_model": "hill_free_n"}))
    out_args, paths = cli_outputs(tmp_path, ".parquet")
    assert main(["--params", str(params_path), "--reads", str(reads_path), *out_args]) == 0
    pc = pl.read_parquet(paths["per_clonotype"])
    h1 = pc.filter(pl.col("clonotypeKey") == "H1")
    assert h1["affinityClass"][0] == "Failed"
    assert h1["fitFailureReason"][0] == "non_monotonic_signal"
    # R17: null-Kd Failed rows must surface at the plot sentinels so the scatter doesn't drop them.
    assert h1["kd"][0] is None
    assert h1["kdPlotPosition"][0] == max(concs) * 10.0
    assert h1["hillPlotPosition"][0] == -1.0


def test_cli_fit_summary_counts_failed_clonotypes(tmp_path):
    # The block shows "All clonotypes failed" from this file; a non-binder is Failed (above range).
    rows = _cli_reads("BINDER", 1e-8, _DEFAULT_MODEL_CONCS, [1, 2, 3, 4])
    rows += _cli_reads("FLAT", 1e2, _DEFAULT_MODEL_CONCS, [1, 2, 3, 4])
    reads_path = tmp_path / "reads.tsv"
    pl.DataFrame(rows).write_csv(reads_path, separator="\t")
    out_args, paths = cli_outputs(tmp_path)
    summary_path = tmp_path / "fit_summary.json"
    assert main(["--reads", str(reads_path), *out_args, *_CELLS_ARGS, "--out-fit-summary", str(summary_path)]) == 0
    pc = pl.read_csv(paths["per_clonotype"], separator="\t")
    expected_failed = pc.filter(pl.col("affinityClass") == "Failed").height
    assert json.loads(summary_path.read_text()) == {"clonotypes": 2, "failed": expected_failed}
    assert expected_failed >= 1
