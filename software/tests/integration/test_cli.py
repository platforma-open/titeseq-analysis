"""CLI smoke tests for main.main(argv) — exercise argparse, file I/O, output writers."""

from __future__ import annotations

import json
import math

import numpy as np
import polars as pl

from main import main


def _hill_reads(clonotype: str, true_kd: float, concs, bins, per_conc=500):
    baseline = 1.5
    top = baseline + math.exp(math.log(2.0))
    rows = []
    for i, c in enumerate(concs):
        target = baseline + (top - baseline) * c / (true_kd + c)
        weights = np.exp(-0.5 * ((np.array(bins, dtype=float) - target) / 0.35) ** 2)
        weights /= weights.sum()
        counts = np.round(weights * per_conc).astype(int)
        if counts.sum() == 0:
            counts[len(counts) // 2] = 1
        for j, b in enumerate(bins):
            rows.append(
                {
                    "clonotypeKey": clonotype,
                    "sampleId": f"s_c{i}_b{b}",
                    "concentrationStr": str(c),
                    "concentration": float(c),
                    "bin": int(b),
                    "reads": int(counts[j]),
                }
            )
    return rows


# Parquet round-trip — verifies argparse wiring, file reader, and parquet writer branch
def test_cli_parquet_roundtrip(tmp_path):
    rows = _hill_reads("G1", 1e-8, [1e-10, 1e-9, 1e-8, 1e-7, 1e-6], [1, 2, 3, 4])
    reads_path = tmp_path / "reads.parquet"
    pl.DataFrame(rows).write_parquet(reads_path)

    pc_path = tmp_path / "pc.parquet"
    mb_path = tmp_path / "mb.parquet"
    fmb_path = tmp_path / "fmb.parquet"

    rc = main(
        [
            "--reads",
            str(reads_path),
            "--out-per-clonotype",
            str(pc_path),
            "--out-mean-bin",
            str(mb_path),
            "--out-fitted-mean-bin",
            str(fmb_path),
        ]
    )
    assert rc == 0
    assert pc_path.exists() and mb_path.exists() and fmb_path.exists()

    pc = pl.read_parquet(pc_path)
    assert pc.height == 1
    assert pc["clonotypeKey"][0] == "G1"


# TSV branch — exercises the non-parquet writer path
def test_cli_tsv_output(tmp_path):
    rows = _hill_reads("G1", 1e-8, [1e-10, 1e-9, 1e-8, 1e-7, 1e-6], [1, 2, 3, 4])
    reads_path = tmp_path / "reads.parquet"
    pl.DataFrame(rows).write_parquet(reads_path)

    pc_path = tmp_path / "pc.tsv"
    mb_path = tmp_path / "mb.tsv"
    fmb_path = tmp_path / "fmb.tsv"

    rc = main(
        [
            "--reads",
            str(reads_path),
            "--out-per-clonotype",
            str(pc_path),
            "--out-mean-bin",
            str(mb_path),
            "--out-fitted-mean-bin",
            str(fmb_path),
        ]
    )
    assert rc == 0
    assert pc_path.read_text().startswith("clonotypeKey")


# --params JSON — ensures FitParams override path is wired through argparse
def test_cli_with_params_json(tmp_path):
    rows = _hill_reads("G1", 1e-8, [1e-10, 1e-9, 1e-8, 1e-7, 1e-6], [1, 2, 3, 4])
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
            }
        )
    )

    rc = main(
        [
            "--reads",
            str(reads_path),
            "--out-per-clonotype",
            str(tmp_path / "pc.parquet"),
            "--out-mean-bin",
            str(tmp_path / "mb.parquet"),
            "--out-fitted-mean-bin",
            str(tmp_path / "fmb.parquet"),
            "--params",
            str(params_path),
        ]
    )
    assert rc == 0


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

    pc_path = tmp_path / "pc.parquet"
    rc = main(
        [
            "--reads",
            str(reads_path),
            "--out-per-clonotype",
            str(pc_path),
            "--out-mean-bin",
            str(tmp_path / "mb.parquet"),
            "--out-fitted-mean-bin",
            str(tmp_path / "fmb.parquet"),
        ]
    )
    assert rc == 0
    pc = pl.read_parquet(pc_path)
    h1 = pc.filter(pl.col("clonotypeKey") == "H1")
    assert h1["affinityClass"][0] == "Failed"
    assert h1["fitFailureReason"][0] == "non_monotonic_signal"
    # R17: null-K_D Failed rows must surface at the plot sentinels so the scatter doesn't drop them.
    assert h1["kd"][0] is None
    assert h1["kdPlotPosition"][0] == max(concs) * 10.0
    assert h1["hillPlotPosition"][0] == -1.0
