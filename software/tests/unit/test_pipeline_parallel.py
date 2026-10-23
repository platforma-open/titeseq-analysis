"""Behavioral tests for the parallel vs. serial fit executor in pipeline._execute_fits.

The pool is a performance detail, but the contract is strict: parallel and serial
paths must produce identical results for the same inputs. These tests exercise
both paths and compare outputs.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from functools import partial

import numpy as np
import pytest
from helpers import hill_curve

import pipeline
from hill_fit import FitResult, fit_one_clonotype


def _make_tasks(n_tasks: int, rng: np.random.Generator) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    """Generate n_tasks realistic (x, y, w) triples for the Hill fitter."""
    concs = np.array([0.0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6])
    xs, ys, ws = [], [], []
    for _ in range(n_tasks):
        # Spread Kd across the grid to keep fits non-degenerate.
        kd = 10 ** rng.uniform(-10, -7)
        baseline = 1.0
        amplitude = math.log(rng.uniform(1.0, 3.0))
        y = hill_curve(concs, baseline, amplitude, kd, n=1.0)
        # Small noise, well under δ_bin — doesn't push survivors below the gate.
        y = y + rng.normal(0, 0.02, size=y.shape)
        xs.append(concs.copy())
        ys.append(y)
        ws.append(np.ones_like(concs) * 100.0)
    return xs, ys, ws


def _worker() -> Callable[[np.ndarray, np.ndarray, np.ndarray], FitResult]:
    return partial(
        fit_one_clonotype,
        baseline_fixed=1.0,
        bin_mode=True,
        max_bin_label=8,
    )


class TestExecuteFitsParity:
    """Serial and parallel paths must produce identical FitResults."""

    # 80 survivors: comfortably above the default _PARALLEL_FIT_MIN_SURVIVORS gate (50)
    # so the parallel branch runs, and divides across 4 workers cleanly.
    def test_parallel_matches_serial(self, monkeypatch):
        rng = np.random.default_rng(seed=42)
        xs, ys, ws = _make_tasks(80, rng)
        worker = _worker()

        # Force serial: threshold above n_tasks.
        monkeypatch.setattr(pipeline, "_PARALLEL_FIT_MIN_SURVIVORS", 10_000)
        serial_fits = pipeline._execute_fits(worker, xs, ys, ws, n_survivors=80, step=20)

        # Force parallel: threshold low, and ensure available_cpus > 1 for the branch to trigger.
        monkeypatch.setattr(pipeline, "_PARALLEL_FIT_MIN_SURVIVORS", 1)
        monkeypatch.setattr(pipeline, "available_cpus", lambda: 4)
        parallel_fits = pipeline._execute_fits(worker, xs, ys, ws, n_survivors=80, step=20)

        assert len(serial_fits) == len(parallel_fits) == 80
        for i, (s, p) in enumerate(zip(serial_fits, parallel_fits)):
            # Converged-ness must match; numeric fields within FP tolerance.
            assert s.converged == p.converged, f"task {i}: converged mismatch"
            assert s.reason == p.reason, f"task {i}: reason mismatch"
            if s.converged:
                assert s.kd == pytest.approx(p.kd, rel=1e-9), f"task {i}: kd mismatch"
                assert s.n == pytest.approx(p.n, rel=1e-9), f"task {i}: n mismatch"
                assert s.r2_w == pytest.approx(p.r2_w, rel=1e-9), f"task {i}: r² mismatch"

    # Below the threshold, the helper must stay serial regardless of CPU count
    # (pool startup outweighs savings for tiny batches).
    def test_small_batch_stays_serial(self, monkeypatch):
        rng = np.random.default_rng(seed=7)
        xs, ys, ws = _make_tasks(5, rng)
        worker = _worker()

        # Default threshold (50) with 5 tasks — parallel branch must not activate.
        # We detect the serial path by patching ProcessPoolExecutor to explode if used.
        class _ShouldNotRun:
            def __init__(self, *a: object, **k: object) -> None:
                raise AssertionError("ProcessPoolExecutor must not be used below threshold")

        monkeypatch.setattr(pipeline, "ProcessPoolExecutor", _ShouldNotRun)
        fits = pipeline._execute_fits(worker, xs, ys, ws, n_survivors=5, step=1)
        assert len(fits) == 5

    def test_empty_input_returns_empty(self):
        worker = _worker()
        assert pipeline._execute_fits(worker, [], [], [], n_survivors=0, step=1) == []
