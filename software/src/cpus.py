"""CPUs this process may use: the container's CPU quota, not the host's core count.

`os.cpu_count()` reports every core on the node. Under a cgroup CPU limit (the workflow's
`cpu(4)`) that oversubscribes the grant, and each extra fit worker adds its own memory.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

CGROUP_ROOT = Path("/sys/fs/cgroup")


def _cgroup_quota(root: Path) -> float | None:
    try:
        quota, period = (root / "cpu.max").read_text().split()[:2]  # cgroup v2: "<quota> <period>" or "max <period>"
        return None if quota == "max" else int(quota) / int(period)
    except (OSError, ValueError):
        pass
    try:
        quota_us = int((root / "cpu" / "cpu.cfs_quota_us").read_text())  # cgroup v1; -1 = unlimited
        period_us = int((root / "cpu" / "cpu.cfs_period_us").read_text())
        return None if quota_us <= 0 else quota_us / period_us
    except (OSError, ValueError):
        return None


def _host_cpus() -> int:
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def available_cpus(cgroup_root: Path = CGROUP_ROOT) -> int:
    """At least 1; at most the host's usable cores and the floor of the cgroup CPU quota."""
    n = _host_cpus()
    quota = _cgroup_quota(cgroup_root)
    if quota is not None:
        n = min(n, max(1, math.floor(quota)))
    return max(1, n)
