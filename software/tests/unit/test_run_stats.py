from __future__ import annotations

import math

import polars as pl
import pytest

from run_stats import fit_reads, run_stats, with_fit_reads

BOUNDS = (1e-10, 1e-6)


def _frame(rows):
    return pl.DataFrame(
        rows,
        schema={"kd": pl.Float64, "affinityClass": pl.Utf8, "kdOutOfRange": pl.Boolean, "fitReads": pl.Int64},
        orient="row",
    )


def test_fit_reads_sums_each_key():
    reads = pl.DataFrame({"k": ["a", "a", "b"], "reads": [3, None, 5]})
    assert dict(fit_reads(reads, "k", "reads").iter_rows()) == {"a": 3, "b": 5}


def test_a_key_with_no_reads_gets_zero():
    frame = with_fit_reads(pl.DataFrame({"k": ["a", "b"]}), "k", {"a": 7})
    assert frame["fitReads"].to_list() == [7, 0]


def test_classes_walls_and_median_reads():
    frame = _frame(
        [
            (1e-8, "Good", False, 100),
            (1e-8, "Good", False, 300),
            (1e-10, "Below range", True, 50),
            (1e-6, "Above range", True, 20),
            (None, "Failed", None, 2),
        ]
    )
    stats = run_stats(frame, grain="variant", conc_bounds=BOUNDS, concentrations=9, replicates=["r1", "r2"])
    assert stats["classes"] == {"Good": 2, "Below range": 1, "Above range": 1, "Failed": 1}
    assert stats["walls"] == {"below": 1, "above": 1}
    assert stats["medianReadsByClass"]["Good"] == 200
    assert stats["concentrations"] == {"count": 9, "log10Lowest": -10.0, "log10Highest": -6.0}
    assert stats["replicates"] == ["r1", "r2"]


def test_the_parent_headroom_is_its_distance_to_each_wall():
    stats = run_stats(
        _frame([]),
        grain="aa",
        conc_bounds=BOUNDS,
        concentrations=9,
        replicates=[],
        parents=[{"parentId": "P", "log10Kd": -8.8, "log10KdSe": 0.1}],
    )
    [p] = stats["parents"]
    assert p["tightRoom"] == pytest.approx(1.2)
    assert p["weakRoom"] == pytest.approx(2.8)
    assert math.isclose(p["log10Kd"], -8.8)
