"""End-to-end bin-mode pipeline tests via main.run()."""

from __future__ import annotations

import math

import numpy as np
import polars as pl
import pytest

from constants import DEFAULT_PARAMS, FitParams
from pipeline import run


def _build_bin_reads_for_hill(
    clonotype: str,
    true_kd: float,
    true_n: float,
    baseline: float,
    amplitude: float,
    concs: list[float],
    conc_strs: list[str] | None,
    bins: list[int],
    reads_per_conc: int,
) -> list[dict]:
    """Generate reads so the computed mean_bin ≈ Hill(c). Deterministic (no noise).

    For each concentration c:
      target_mean_bin = Hill(c)
      Distribute reads_per_conc across bins using a Gaussian centered at target_mean_bin.
    """
    conc_strs = conc_strs or [str(c) for c in concs]
    rows = []
    for i, c in enumerate(concs):
        top = baseline + math.exp(amplitude)
        if c == 0:
            target = baseline
        else:
            target = baseline + (top - baseline) * (c**true_n) / (true_kd**true_n + c**true_n)
        centers = np.array(bins, dtype=float)
        weights = np.exp(-0.5 * ((centers - target) / 0.35) ** 2)
        weights /= weights.sum()
        counts = np.round(weights * reads_per_conc).astype(int)
        # Ensure at least 1 read total.
        if counts.sum() == 0:
            counts[len(counts) // 2] = 1
        for j, b in enumerate(bins):
            rows.append(
                {
                    "clonotypeKey": clonotype, "sampleId": f"s_c{i}_b{b}",
                    "concentrationStr": conc_strs[i], "concentration": float(c),
                    "bin": int(b), "reads": int(counts[j]),
                }
            )
    return rows


class TestBinModePipeline:
    # Noiseless Hill reads → Good class, K_D recovered within relative tolerance.
    def test_good_clonotype_recovered(self):
        rows = _build_bin_reads_for_hill(
            "G1", true_kd=10.0, true_n=1.0, baseline=1.5,
            amplitude=math.log(2.0),
            concs=[0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0],
            conc_strs=None,
            bins=[1, 2, 3, 4], reads_per_conc=500,
        )
        reads = pl.DataFrame(rows)
        out = run(reads, params=DEFAULT_PARAMS)
        pc = out["per_clonotype"].filter(pl.col("clonotypeKey") == "G1")
        assert pc["affinityClass"][0] in {"Good", "Partial"}
        # Synthesized-bin quantization distorts K_D; check order of magnitude only.
        assert pc["kd"][0] == pytest.approx(10.0, rel=2.0)

    # c=0 control anchors baseline → 3-param fit path exercised (R6).
    def test_c0_control_fixes_baseline(self):
        concs = [0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
        rows = _build_bin_reads_for_hill(
            "G1", true_kd=5.0, true_n=1.0, baseline=1.0,
            amplitude=math.log(2.5),
            concs=concs, conc_strs=None,
            bins=[1, 2, 3, 4], reads_per_conc=500,
        )
        reads = pl.DataFrame(rows)
        out = run(reads, params=DEFAULT_PARAMS)
        pc = out["per_clonotype"].filter(pl.col("clonotypeKey") == "G1")
        assert pc["affinityClass"][0] in {"Good", "Partial"}
        assert pc["kd"][0] == pytest.approx(5.0, rel=2.0)

    # Output table shape: fittedMeanBin only emitted for converged fits.
    def test_fitted_mean_bin_only_for_converged(self):
        rows = _build_bin_reads_for_hill(
            "G1", true_kd=10.0, true_n=1.0, baseline=1.5,
            amplitude=math.log(2.0),
            concs=[0.1, 1.0, 10.0, 100.0, 1000.0],
            conc_strs=None, bins=[1, 2, 3, 4], reads_per_conc=500,
        )
        rows += _build_bin_reads_for_hill(
            "F1", true_kd=0.0, true_n=1.0, baseline=2.0, amplitude=math.log(0.001),
            concs=[0.1, 1.0, 10.0, 100.0, 1000.0],
            conc_strs=None, bins=[1, 2, 3, 4], reads_per_conc=500,
        )
        reads = pl.DataFrame(rows)
        out = run(reads)
        pc = out["per_clonotype"]
        fmb = out["fitted_mean_bin"]
        # F1 has amplitude ≈ log(0.001) → δ-fail (convergence_failure)
        assert pc.filter(pl.col("clonotypeKey") == "F1")["affinityClass"][0] == "Failed"
        assert fmb.filter(pl.col("clonotypeKey") == "F1").height == 0
        assert fmb.filter(pl.col("clonotypeKey") == "G1").height > 0

    # R14: meanBin output excludes c=0 rows.
    def test_mean_bin_output_excludes_c0(self):
        concs = [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0]
        rows = _build_bin_reads_for_hill(
            "G1", true_kd=10.0, true_n=1.0, baseline=1.0,
            amplitude=math.log(2.0),
            concs=concs, conc_strs=None, bins=[1, 2, 3, 4], reads_per_conc=400,
        )
        reads = pl.DataFrame(rows)
        out = run(reads)
        mb = out["mean_bin"].filter(pl.col("clonotypeKey") == "G1")
        assert "0" not in mb["concentrationStr"].to_list()


class TestAllFailedDataset:
    # Spec edge case: every clonotype fails → block completes without error.
    def test_all_clonotypes_failed_completes(self):
        # All flat signal (below δ) → all convergence_failure.
        concs = [0.1, 1.0, 10.0, 100.0, 1000.0]
        rows = []
        for k in range(3):
            clonotype = f"C{k}"
            for i, c in enumerate(concs):
                for b in [1, 2, 3, 4]:
                    rows.append({
                        "clonotypeKey": clonotype, "sampleId": f"s_c{i}_b{b}",
                        "concentrationStr": str(c), "concentration": float(c),
                        "bin": b, "reads": 25,  # equal across bins → mean_bin = 2.5 flat
                    })
        reads = pl.DataFrame(rows)
        out = run(reads)
        pc = out["per_clonotype"]
        assert pc.height == 3
        assert (pc["affinityClass"] == "Failed").all()
        # Spec requires output file still written; kd column all null.
        assert pc["kd"].null_count() == 3


class TestValidatorWarnings:
    # R2: c=0 without a bin is ambiguous — validator emits a warning to stderr
    # (a non-fatal signal; the row is still processed).
    def test_c0_without_bin_emits_stderr_warning(self, capsys):
        concs = [0.0, 0.1, 1.0, 10.0, 100.0]
        rows = _build_bin_reads_for_hill(
            "G1", true_kd=10.0, true_n=1.0, baseline=1.0,
            amplitude=math.log(2.0),
            concs=concs, conc_strs=None, bins=[1, 2, 3, 4], reads_per_conc=400,
        )
        # Inject a c=0 row with bin=None (ambiguous control)
        rows.append({
            "clonotypeKey": "G1", "sampleId": "s_c0_unbinned",
            "concentrationStr": "0", "concentration": 0.0,
            "bin": None, "reads": 100,
        })
        reads = pl.DataFrame(rows, schema_overrides={"bin": pl.Int64})
        run(reads)
        assert "ambiguous" in capsys.readouterr().err


class TestInsufficientReads:
    # Low-read clonotype hits insufficient_reads.
    def test_all_below_floor_insufficient_reads(self):
        concs = [0.1, 1.0, 10.0, 100.0, 1000.0]
        rows = []
        for i, c in enumerate(concs):
            for b in [1, 2, 3, 4]:
                rows.append({
                    "clonotypeKey": "Low", "sampleId": f"s_c{i}_b{b}",
                    "concentrationStr": str(c), "concentration": float(c),
                    "bin": b, "reads": 0,  # no reads at all
                })
                # Filler reads under a second clonotype to create non-zero depth
                rows.append({
                    "clonotypeKey": "Filler", "sampleId": f"s_c{i}_b{b}",
                    "concentrationStr": str(c), "concentration": float(c),
                    "bin": b, "reads": 100,
                })
        reads = pl.DataFrame(rows)
        out = run(reads)
        low = out["per_clonotype"].filter(pl.col("clonotypeKey") == "Low")
        assert low["affinityClass"][0] == "Failed"
        assert low["fitFailureReason"][0] == "insufficient_reads"
