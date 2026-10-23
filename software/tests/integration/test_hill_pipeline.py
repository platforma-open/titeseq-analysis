"""End-to-end tests of the opt-in linear Hill model (hill_free_n) through pipeline.run."""

from __future__ import annotations

import math

import polars as pl
from helpers import LINEAR_HILL, hill_reads

from pipeline import run

# Sub-µM grid (0.1 nM to 1 µM), matching realistic Tite-Seq dose ranges.
CONCS = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
BINS = [1, 2, 3, 4]


# Spec edge case: every clonotype fails → block completes without error.
def test_all_clonotypes_failed_completes():
    # Equal reads in every bin give a flat mean bin of 2.5 (below δ) → all convergence_failure.
    rows = [
        {
            "clonotypeKey": f"C{k}",
            "sampleId": f"s_c{i}_b{b}",
            "concentrationStr": str(c),
            "concentration": float(c),
            "bin": b,
            "reads": 25,
        }
        for k in range(3)
        for i, c in enumerate(CONCS)
        for b in BINS
    ]
    pc = run(pl.DataFrame(rows), params=LINEAR_HILL)["per_clonotype"]
    assert pc.height == 3
    assert (pc["affinityClass"] == "Failed").all()
    # Spec requires output file still written; kd column all null.
    assert pc["kd"].null_count() == 3


# FACS sort-fraction correction, end to end. Bin 1 is over-represented at low concentrations and bin 4
# at high ones, the pattern the Adams-Mora-Walczak-Kinney correction targets. The run must keep the
# output schema, still fit, and move mean bin away from the uncorrected run (else it is a no-op).
def test_pipeline_with_sort_fraction():
    rows = hill_reads("G1", 1e-8, CONCS, BINS, baseline=1.0, amplitude=math.log(2.5), per_conc=800)
    sort_fractions_by_conc = {
        1e-10: [0.50, 0.25, 0.15, 0.10],
        1e-9: [0.40, 0.30, 0.20, 0.10],
        1e-8: [0.25, 0.25, 0.25, 0.25],
        1e-7: [0.10, 0.20, 0.30, 0.40],
        1e-6: [0.10, 0.15, 0.25, 0.50],
    }
    reads_plain = pl.DataFrame(rows)
    reads_with_sort = pl.DataFrame(
        [r | {"sort_fraction": sort_fractions_by_conc[r["concentration"]][r["bin"] - 1]} for r in rows]
    )

    out_uncorrected = run(reads_plain, params=LINEAR_HILL)
    out_corrected = run(reads_with_sort, params=LINEAR_HILL, sort_fraction_col="sort_fraction")

    # The correction stays on the existing meanBin column: no new columns or rows.
    assert set(out_corrected["mean_bin"].columns) == set(out_uncorrected["mean_bin"].columns)
    assert out_corrected["mean_bin"].height == out_uncorrected["mean_bin"].height
    assert set(out_corrected["per_clonotype"].columns) == set(out_uncorrected["per_clonotype"].columns)
    assert out_corrected["per_clonotype"].height == 1
    assert not out_corrected["per_clonotype"]["kd"].is_null().all()
    corrected = out_corrected["mean_bin"].sort("concentrationStr")["meanBin"].to_list()
    uncorrected = out_uncorrected["mean_bin"].sort("concentrationStr")["meanBin"].to_list()
    diffs = [abs(c - u) for c, u in zip(corrected, uncorrected) if c is not None and u is not None]
    assert any(d > 0.05 for d in diffs), "sort-fraction correction left mean_bin unchanged — is the kwarg wired?"


# R2: c=0 without a bin is ambiguous — validator emits a WARN to stdout so the
# Tengo workflow's saveStdoutStream() surfaces it in the Fit Log UI.
def test_c0_without_bin_emits_stdout_warning(capsys):
    rows = hill_reads("G1", 1e-8, [0.0, *CONCS[:4]], BINS, baseline=1.0, per_conc=400)
    rows.append(
        {
            "clonotypeKey": "G1",
            "sampleId": "s_c0_unbinned",
            "concentrationStr": "0",
            "concentration": 0.0,
            "bin": None,
            "reads": 100,
        }
    )
    run(pl.DataFrame(rows, schema_overrides={"bin": pl.Int64}), params=LINEAR_HILL)
    captured = capsys.readouterr()
    assert "ambiguous" in captured.out
    assert "ambiguous" not in captured.err


# R4: the same clonotype measured against two antigens is fitted on the target antigen's reads only.
def test_antigen_filter_keeps_only_target():
    rows = hill_reads("C1", 5e-9, CONCS, BINS, antigen="X") + hill_reads("C1", 5e-7, CONCS, BINS, antigen="Y")
    out = run(pl.DataFrame(rows), params=LINEAR_HILL, target_antigen="X", antigen_column_ref="antigen")
    pc = out["per_clonotype"].filter(pl.col("clonotypeKey") == "C1")
    # With only antigen X present, Kd ≈ 5 nM (not 500 nM)
    assert pc["kd"][0] is not None
    assert pc["kd"][0] < 5e-8
