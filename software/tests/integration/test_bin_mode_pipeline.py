"""End-to-end bin-mode pipeline tests via main.run()."""

from __future__ import annotations

import math

import numpy as np
import polars as pl

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
                    "clonotypeKey": clonotype,
                    "sampleId": f"s_c{i}_b{b}",
                    "concentrationStr": conc_strs[i],
                    "concentration": float(c),
                    "bin": int(b),
                    "reads": int(counts[j]),
                }
            )
    return rows


class TestAllFailedDataset:
    # Spec edge case: every clonotype fails → block completes without error.
    def test_all_clonotypes_failed_completes(self):
        # All flat signal (below δ) → all convergence_failure.
        # Sub-µM range (0.1 nM → 1 µM) matching realistic TiteSeq dose grids.
        concs = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
        rows = []
        for k in range(3):
            clonotype = f"C{k}"
            for i, c in enumerate(concs):
                for b in [1, 2, 3, 4]:
                    rows.append(
                        {
                            "clonotypeKey": clonotype,
                            "sampleId": f"s_c{i}_b{b}",
                            "concentrationStr": str(c),
                            "concentration": float(c),
                            "bin": b,
                            "reads": 25,  # equal across bins → mean_bin = 2.5 flat
                        }
                    )
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
        # Sub-µM grid; K_D placed in the middle of the dose range (10 nM).
        concs = [0.0, 1e-10, 1e-9, 1e-8, 1e-7]
        rows = _build_bin_reads_for_hill(
            "G1",
            true_kd=1e-8,
            true_n=1.0,
            baseline=1.0,
            amplitude=math.log(2.0),
            concs=concs,
            conc_strs=None,
            bins=[1, 2, 3, 4],
            reads_per_conc=400,
        )
        # Inject a c=0 row with bin=None (ambiguous control)
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
        reads = pl.DataFrame(rows, schema_overrides={"bin": pl.Int64})
        run(reads)
        assert "ambiguous" in capsys.readouterr().err


