"""available_cpus: fit worker pools must follow the container's CPU quota, not the node's cores."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from cpus import _host_cpus, available_cpus


def _v2(root: Path, text: str) -> Path:
    (root / "cpu.max").write_text(text)
    return root


def _v1(root: Path, quota_us: int, period_us: int) -> Path:
    (root / "cpu").mkdir()
    (root / "cpu" / "cpu.cfs_quota_us").write_text(f"{quota_us}\n")
    (root / "cpu" / "cpu.cfs_period_us").write_text(f"{period_us}\n")
    return root


@pytest.mark.parametrize(
    ("cpu_max", "expected_cap"),
    [
        ("400000 100000\n", 4),  # the workflow's cpu(4)
        ("150000 100000\n", 1),  # fractional quota rounds down
        ("50000 100000\n", 1),  # below one CPU still runs one worker
    ],
)
def test_cgroup_v2_quota_caps_workers(tmp_path: Path, cpu_max: str, expected_cap: int) -> None:
    assert available_cpus(_v2(tmp_path, cpu_max)) == min(expected_cap, _host_cpus())


def test_cgroup_v2_unlimited_uses_host_cores(tmp_path: Path) -> None:
    assert available_cpus(_v2(tmp_path, "max 100000\n")) == _host_cpus()


def test_cgroup_v1_quota_caps_workers(tmp_path: Path) -> None:
    assert available_cpus(_v1(tmp_path, 200000, 100000)) == min(2, _host_cpus())


def test_cgroup_v1_unlimited_uses_host_cores(tmp_path: Path) -> None:
    assert available_cpus(_v1(tmp_path, -1, 100000)) == _host_cpus()


def test_no_cgroup_files_uses_host_cores(tmp_path: Path) -> None:
    # Local backend on a laptop: no cgroup CPU controller is visible.
    assert available_cpus(tmp_path) == _host_cpus()


def test_unparseable_cpu_max_is_ignored(tmp_path: Path) -> None:
    assert available_cpus(_v2(tmp_path, "garbage\n")) == _host_cpus()


@given(quota=st.integers(min_value=1, max_value=10**7), period=st.integers(min_value=1000, max_value=10**6))
def test_workers_never_exceed_quota_or_host(tmp_path_factory: pytest.TempPathFactory, quota: int, period: int) -> None:
    root = _v2(tmp_path_factory.mktemp("cg"), f"{quota} {period}\n")
    n = available_cpus(root)
    assert 1 <= n <= _host_cpus()
    assert n <= max(1, math.floor(quota / period))
