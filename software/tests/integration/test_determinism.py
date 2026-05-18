"""Identical inputs must produce byte-identical outputs across runs."""

from __future__ import annotations

import filecmp
import math
from pathlib import Path

import numpy as np
import polars as pl

from main import main


def _hill_reads(clonotype: str, true_kd: float, concs, bins, per_conc=500):
    """Synthetic per-clonotype Hill-binding reads — mirrors test_cli helper."""
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


def _run_cli(reads_path: Path, out_dir: Path) -> dict[str, Path]:
    pc = out_dir / "per_clonotype.tsv"
    mb = out_dir / "mean_bin.tsv"
    fmb = out_dir / "fitted_mean_bin.tsv"
    cv = out_dir / "concentration_value.tsv"
    rc = main(
        [
            "--reads", str(reads_path),
            "--out-per-clonotype", str(pc),
            "--out-mean-bin", str(mb),
            "--out-fitted-mean-bin", str(fmb),
            "--out-concentration-value", str(cv),
        ]
    )
    assert rc == 0
    return {"per_clonotype": pc, "mean_bin": mb, "fitted_mean_bin": fmb, "concentration_value": cv}


def test_outputs_are_byte_identical_across_runs(tmp_path):
    """Run main() twice on the same fixture, byte-compare all four output TSVs.

    Multiple clonotypes and multiple concentrations make the test sensitive to
    non-deterministic polars.group_by, dict iteration order, and any other
    ordering hazard that could surface as a CIDConflict downstream.
    """
    concs = [1e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 1e-6]
    bins = [1, 2, 3, 4]
    rows: list[dict] = []
    for clonotype, kd in [("CLONE_A", 1e-8), ("CLONE_B", 3e-9), ("CLONE_C", 5e-9), ("CLONE_D", 1e-7)]:
        rows.extend(_hill_reads(clonotype, kd, concs, bins))

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
