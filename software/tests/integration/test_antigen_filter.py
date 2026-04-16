"""R4: antigen filtering end-to-end (target_antigen is kept, others discarded)."""

from __future__ import annotations

import math

import numpy as np
import polars as pl

from main import run


def _good_reads_for(clonotype, antigen, true_kd, baseline, concs, bins, per_conc=500):
    rows = []
    top = baseline + math.exp(math.log(2.0))
    for i, c in enumerate(concs):
        target = baseline + (top - baseline) * c / (true_kd + c)
        weights = np.exp(-0.5 * ((np.array(bins, dtype=float) - target) / 0.35) ** 2)
        weights /= weights.sum()
        counts = np.round(weights * per_conc).astype(int)
        for j, b in enumerate(bins):
            rows.append({
                "clonotypeKey": clonotype, "sampleId": f"{antigen}_s_c{i}_b{b}",
                "concentrationStr": str(c), "concentration": float(c), "bin": int(b),
                "reads": int(counts[j]), "antigen": antigen,
            })
    return rows


def test_antigen_filter_keeps_only_target():
    concs = [0.1, 1.0, 10.0, 100.0, 1000.0]
    bins = [1, 2, 3, 4]
    rows = _good_reads_for("C1", "X", true_kd=5.0, baseline=1.5, concs=concs, bins=bins)
    # Same clonotype under a different antigen with a very different K_D
    rows += _good_reads_for("C1", "Y", true_kd=500.0, baseline=1.5, concs=concs, bins=bins)
    reads = pl.DataFrame(rows)

    out_x = run(reads, target_antigen="X", antigen_column_ref="antigen")
    pc = out_x["per_clonotype"].filter(pl.col("clonotypeKey") == "C1")
    # With only antigen X present, K_D ≈ 5 (not 500)
    assert pc["kd"][0] is not None
    assert pc["kd"][0] < 50.0
