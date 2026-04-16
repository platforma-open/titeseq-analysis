"""Parallel fit path must produce identical output to the serial path."""

from __future__ import annotations

import math

import numpy as np
import polars as pl

import pipeline
from pipeline import run


def _noiseless_bin_reads(clonotype: str, kd: float, n: float) -> list[dict]:
    baseline, amplitude = 1.0, math.log(3.0)
    concs = [0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0]
    bins = [1, 2, 3, 4]
    rows: list[dict] = []
    top = baseline + math.exp(amplitude)
    for i, c in enumerate(concs):
        if c == 0:
            target = baseline
        else:
            target = baseline + (top - baseline) * (c**n) / (kd**n + c**n)
        centers = np.array(bins, dtype=float)
        weights = np.exp(-0.5 * ((centers - target) / 0.35) ** 2)
        weights /= weights.sum()
        counts = np.round(weights * 100).astype(int)
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


# Serial and parallel must be numerically identical — the pool only changes
# scheduling, not arithmetic. Guards against accidental pickle-related state
# drift in fit_one_clonotype (e.g., hidden globals, RNG seeds).
def test_parallel_matches_serial(monkeypatch):
    # Enough clonotypes to be worth parallelising; threshold lowered so the pool
    # branch actually runs on this fixture.
    monkeypatch.setattr(pipeline, "_PARALLEL_MIN_TASKS", 4)

    kd_grid = [0.5, 2.0, 10.0, 50.0, 200.0]
    rows: list[dict] = []
    for i, kd in enumerate(kd_grid):
        rows += _noiseless_bin_reads(f"C{i}", kd=kd, n=1.0)
    reads = pl.DataFrame(rows)

    serial = run(reads, workers=1)
    parallel = run(reads, workers=2)

    for key in ("per_clonotype", "mean_bin", "fitted_mean_bin"):
        s = serial[key].sort(serial[key].columns)
        p = parallel[key].sort(parallel[key].columns)
        assert s.equals(p), f"{key} frame diverged between serial and parallel"


# Small input stays serial even when workers>1 — the pool threshold protects
# against paying spawn cost on tiny datasets.
def test_pool_not_engaged_below_threshold(monkeypatch):
    called = {"n": 0}
    original = pipeline.ProcessPoolExecutor

    class _Sentinel:
        def __init__(self, *a, **kw):
            called["n"] += 1
            raise AssertionError("ProcessPoolExecutor should not be constructed")

    monkeypatch.setattr(pipeline, "ProcessPoolExecutor", _Sentinel)
    monkeypatch.setattr(pipeline, "_PARALLEL_MIN_TASKS", 1000)

    reads = pl.DataFrame(_noiseless_bin_reads("C0", kd=10.0, n=1.0))
    out = run(reads, workers=4)  # below threshold → serial
    assert called["n"] == 0
    assert out["per_clonotype"].height == 1
    _ = original  # silence unused
