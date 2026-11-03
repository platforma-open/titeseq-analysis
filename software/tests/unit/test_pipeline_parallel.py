"""Behavioral tests for the parallel vs. serial fit executor in pipeline._execute_fits.

The pool is a performance detail, but the contract is strict: parallel and serial
paths must produce identical results for the same inputs. These tests exercise
both paths and compare outputs.
"""

from __future__ import annotations

import math
from functools import partial

import numpy as np
import pytest

import pipeline
from hill_fit import fit_one_clonotype


def _hill_truth(x, baseline, amplitude, kd, n):
    top = baseline + math.exp(amplitude)
    return baseline + (top - baseline) * (x**n) / (kd**n + x**n)


def _make_tasks(n_tasks: int, rng: np.random.Generator):
    """Generate n_tasks realistic (x, y, w) triples for the Hill fitter."""
    concs = np.array([0.0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6])
    xs, ys, ws = [], [], []
    for i in range(n_tasks):
        # Spread Kd across the grid to keep fits non-degenerate.
        kd = 10 ** rng.uniform(-10, -7)
        baseline = 1.0
        amplitude = math.log(rng.uniform(1.0, 3.0))
        y = _hill_truth(concs, baseline, amplitude, kd, n=1.0)
        # Small noise, well under δ_bin — doesn't push survivors below the gate.
        y = y + rng.normal(0, 0.02, size=y.shape)
        xs.append(concs.copy())
        ys.append(y)
        ws.append(np.ones_like(concs) * 100.0)
    return xs, ys, ws


def _worker():
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

        # Force parallel: threshold low, and ensure more than one CPU for the branch to trigger.
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
            def __init__(self, *a, **k):
                raise AssertionError("ProcessPoolExecutor must not be used below threshold")

        monkeypatch.setattr(pipeline, "ProcessPoolExecutor", _ShouldNotRun)
        fits = pipeline._execute_fits(worker, xs, ys, ws, n_survivors=5, step=1)
        assert len(fits) == 5

    def test_empty_input_returns_empty(self):
        worker = _worker()
        assert pipeline._execute_fits(worker, [], [], [], n_survivors=0, step=1) == []

    def test_pool_size_follows_granted_cpus_not_host_cores(self, monkeypatch):
        rng = np.random.default_rng(seed=3)
        xs, ys, ws = _make_tasks(60, rng)
        sizes: list[int] = []

        class _RecordingPool(pipeline.ProcessPoolExecutor):
            def __init__(self, max_workers=None, **kwargs):
                sizes.append(max_workers)
                super().__init__(max_workers=max_workers, **kwargs)

        monkeypatch.setattr(pipeline, "ProcessPoolExecutor", _RecordingPool)
        monkeypatch.setattr(pipeline.os, "cpu_count", lambda: 96)
        monkeypatch.setattr(pipeline, "available_cpus", lambda: 4)
        pipeline._execute_fits(_worker(), xs, ys, ws, n_survivors=60, step=20)
        assert sizes == [4]


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


class TestAvailableCpus:
    """The CPUs a job may use: the container's CPU quota, never the host's core count."""

    @pytest.fixture(autouse=True)
    def _host_has_96_cores(self, monkeypatch):
        monkeypatch.setattr(pipeline, "_host_cpus", lambda: 96)

    def test_cgroup_v2_quota(self, tmp_path):
        _write(tmp_path / "cpu.max", "400000 100000\n")
        assert pipeline.available_cpus(tmp_path) == 4

    def test_cgroup_v2_fractional_quota_rounds_up(self, tmp_path):
        _write(tmp_path / "cpu.max", "250000 100000\n")
        assert pipeline.available_cpus(tmp_path) == 3

    def test_cgroup_v2_unlimited_falls_back_to_host(self, tmp_path):
        _write(tmp_path / "cpu.max", "max 100000\n")
        assert pipeline.available_cpus(tmp_path) == 96

    def test_cgroup_v1_quota(self, tmp_path):
        _write(tmp_path / "cpu" / "cpu.cfs_quota_us", "200000\n")
        _write(tmp_path / "cpu" / "cpu.cfs_period_us", "100000\n")
        assert pipeline.available_cpus(tmp_path) == 2

    def test_cgroup_v1_unlimited_falls_back_to_host(self, tmp_path):
        _write(tmp_path / "cpu" / "cpu.cfs_quota_us", "-1\n")
        _write(tmp_path / "cpu" / "cpu.cfs_period_us", "100000\n")
        assert pipeline.available_cpus(tmp_path) == 96

    def test_no_cgroup_files_falls_back_to_host(self, tmp_path):
        assert pipeline.available_cpus(tmp_path) == 96

    @pytest.mark.parametrize("text", ["max\n", "\n", "abc 100000\n", "100000 0\n", "1 2 3\n"])
    def test_a_malformed_cgroup_v2_file_falls_back_to_host(self, tmp_path, text):
        _write(tmp_path / "cpu.max", text)
        assert pipeline.available_cpus(tmp_path) == 96

    @pytest.mark.parametrize("quota,period", [("abc\n", "100000\n"), ("200000\n", "0\n")])
    def test_a_malformed_cgroup_v1_file_falls_back_to_host(self, tmp_path, quota, period):
        _write(tmp_path / "cpu" / "cpu.cfs_quota_us", quota)
        _write(tmp_path / "cpu" / "cpu.cfs_period_us", period)
        assert pipeline.available_cpus(tmp_path) == 96

    def test_an_unreadable_cgroup_file_falls_back_to_host(self, tmp_path, monkeypatch):
        _write(tmp_path / "cpu.max", "200000 100000\n")

        def refuse(self, *args, **kwargs):
            raise PermissionError(13, "Permission denied", str(self))

        # Patched rather than chmod-ed: a test run as root reads a mode-000 file anyway.
        monkeypatch.setattr(pipeline.Path, "read_text", refuse)
        assert pipeline.available_cpus(tmp_path) == 96

    def test_quota_above_host_is_capped_at_host(self, tmp_path, monkeypatch):
        monkeypatch.setattr(pipeline, "_host_cpus", lambda: 2)
        _write(tmp_path / "cpu.max", "800000 100000\n")
        assert pipeline.available_cpus(tmp_path) == 2
