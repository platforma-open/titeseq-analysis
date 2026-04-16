"""CLI worker-count parsing (--workers)."""

from __future__ import annotations

import os

import pytest

from main import _resolve_workers


# None means "no flag given" — default to serial (1).
def test_none_defaults_to_one():
    assert _resolve_workers(None) == 1


# 'auto' expands to the host's CPU count; fall back to 1 if cpu_count() returns None.
def test_auto_uses_cpu_count():
    expected = os.cpu_count() or 1
    assert _resolve_workers("auto") == expected


@pytest.mark.parametrize("value, expected", [("1", 1), ("4", 4), ("16", 16)])
def test_int_strings_parse(value, expected):
    assert _resolve_workers(value) == expected


# Non-positive values are user error — raise rather than silently downgrading.
@pytest.mark.parametrize("value", ["0", "-1"])
def test_non_positive_rejected(value):
    with pytest.raises(ValueError, match=">= 1"):
        _resolve_workers(value)
