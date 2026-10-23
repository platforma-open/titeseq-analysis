"""R4: antigen filtering end-to-end (target_antigen is kept, others discarded)."""

from __future__ import annotations

import polars as pl
from helpers import hill_reads

from constants import FitParams
from pipeline import run

# Fixtures are built on the linear Hill model over mean bin, without cell counts or a 0 M control,
# so they exercise the opt-in hill_free_n model.
LINEAR_HILL = FitParams(fit_model="hill_free_n")


def test_antigen_filter_keeps_only_target():
    # Sub-µM grid; Kd values scaled to match (5 nM vs 500 nM antibodies).
    concs = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
    bins = [1, 2, 3, 4]
    rows = hill_reads("C1", 5e-9, concs, bins, antigen="X")
    # Same clonotype under a different antigen with a very different Kd
    rows += hill_reads("C1", 5e-7, concs, bins, antigen="Y")
    reads = pl.DataFrame(rows)

    out_x = run(reads, params=LINEAR_HILL, target_antigen="X", antigen_column_ref="antigen")
    pc = out_x["per_clonotype"].filter(pl.col("clonotypeKey") == "C1")
    # With only antigen X present, Kd ≈ 5 nM (not 500 nM)
    assert pc["kd"][0] is not None
    assert pc["kd"][0] < 5e-8
