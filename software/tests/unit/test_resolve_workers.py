"""CLI worker-count parsing (--workers)."""

from __future__ import annotations

import os

import pytest

from main import _available_cpus, _resolve_workers


# None means "no flag given" — default to serial (1).
def test_none_defaults_to_one():
    assert _resolve_workers(None) == 1


# 'auto' expands to the affinity-aware CPU budget — on Linux this respects
# cgroup/cpuset limits (containers), falling back to os.cpu_count() elsewhere.
def test_auto_uses_available_cpu_budget():
    assert _resolve_workers("auto") == _available_cpus()


# Affinity query, when available, must never exceed os.cpu_count().
# Guards against the helper returning a nonsense value on constrained hosts.
def test_available_cpus_bounded_by_host():
    budget = _available_cpus()
    host = os.cpu_count() or 1
    assert 1 <= budget <= host


@pytest.mark.parametrize("value, expected", [("1", 1), ("4", 4), ("16", 16)])
def test_int_strings_parse(value, expected):
    assert _resolve_workers(value) == expected


# Non-positive values are user error — raise rather than silently downgrading.
@pytest.mark.parametrize("value", ["0", "-1"])
def test_non_positive_rejected(value):
    with pytest.raises(ValueError, match=">= 1"):
        _resolve_workers(value)
