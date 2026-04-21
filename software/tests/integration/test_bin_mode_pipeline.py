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


class TestSortFractionIntegration:
    """End-to-end pipeline test with FACS sort-fraction correction active.

    Uses a mildly-skewed sort where bin 1 is over-represented at low
    concentrations and bin 4 at high concentrations — the same pattern the
    Adams-Mora-Walczak-Kinney correction was designed for. Asserts that the
    pipeline runs to completion, outputs keep their schema, and the Hill fit
    succeeds. A direction-of-shift assertion against the uncorrected run
    guards against the correction being silently a no-op.
    """

    @staticmethod
    def _build_skewed_reads() -> tuple[pl.DataFrame, pl.DataFrame]:
        concs = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
        bins = [1, 2, 3, 4]
        rows = _build_bin_reads_for_hill(
            "G1",
            true_kd=1e-8,
            true_n=1.0,
            baseline=1.0,
            amplitude=math.log(2.5),
            concs=concs,
            conc_strs=None,
            bins=bins,
            reads_per_conc=800,
        )
        # Skew: at the lowest concentration, bin 1 catches 50% of cells; at the
        # highest, bin 4 catches 50%. Intermediate concentrations fall between.
        # Fractions sum to 1.0 per concentration exactly (within float).
        sort_fractions_by_conc: dict[float, list[float]] = {
            1e-10: [0.50, 0.25, 0.15, 0.10],
            1e-9: [0.40, 0.30, 0.20, 0.10],
            1e-8: [0.25, 0.25, 0.25, 0.25],
            1e-7: [0.10, 0.20, 0.30, 0.40],
            1e-6: [0.10, 0.15, 0.25, 0.50],
        }
        with_fraction = []
        for r in rows:
            r2 = dict(r)
            r2["sort_fraction"] = sort_fractions_by_conc[r["concentration"]][r["bin"] - 1]
            with_fraction.append(r2)
        return pl.DataFrame(rows), pl.DataFrame(with_fraction)

    def test_pipeline_with_sort_fraction(self):
        reads_plain, reads_with_sort = self._build_skewed_reads()

        out_uncorrected = run(reads_plain)
        out_corrected = run(reads_with_sort, sort_fraction_col="sort_fraction")

        # Schema must match the uncorrected run — the correction stays on the
        # existing meanBin column, no new columns or rows.
        assert set(out_corrected["mean_bin"].columns) == set(out_uncorrected["mean_bin"].columns)
        assert out_corrected["mean_bin"].height == out_uncorrected["mean_bin"].height
        assert set(out_corrected["per_clonotype"].columns) == set(out_uncorrected["per_clonotype"].columns)

        # Both runs must fit the clonotype — skew doesn't break convergence on
        # this deterministic (no-noise) fixture.
        assert out_corrected["per_clonotype"].height == 1
        assert not out_corrected["per_clonotype"]["kd"].is_null().all()
        # Correction must shift mean_bin away from the uncorrected values at
        # the skewed concentrations — if they're identical we haven't actually
        # wired the correction through.
        corrected_sorted = out_corrected["mean_bin"].sort("concentration")
        uncorrected_sorted = out_uncorrected["mean_bin"].sort("concentration")
        diffs = [
            abs(c - u)
            for c, u in zip(corrected_sorted["meanBin"].to_list(), uncorrected_sorted["meanBin"].to_list())
            if c is not None and u is not None
        ]
        assert any(d > 0.05 for d in diffs), "sort-fraction correction left mean_bin unchanged — is the kwarg wired?"


class TestValidatorWarnings:
    # R2: c=0 without a bin is ambiguous — validator emits a WARN to stdout so the
    # Tengo workflow's saveStdoutStream() surfaces it in the Fit Log UI.
    def test_c0_without_bin_emits_stdout_warning(self, capsys):
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
        captured = capsys.readouterr()
        assert "ambiguous" in captured.out
        assert "ambiguous" not in captured.err


