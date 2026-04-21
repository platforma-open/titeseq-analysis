"""Regression: simulate Poisson reads at the bin level, assert K_D recovery on enough clonotypes.

Spec R-test: synthesize many clonotypes with known K_D, run through the pipeline,
verify >= 80% recovered within a factor (relaxed below from plan's 10% on bin-level reads).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import polars as pl
import pytest

from pipeline import run


@dataclass
class PoissonConfig:
    true_kd: float
    true_n: float
    baseline: float
    amplitude: float
    concs: list[float]
    n_clonotypes: int
    n_bins: int
    reads_per_conc: int
    seed: int


def _hill(x, baseline, amplitude, kd, n):
    top = baseline + math.exp(amplitude)
    return baseline + (top - baseline) * (x**n) / (kd**n + x**n)


def build_poisson_reads(cfg: PoissonConfig) -> pl.DataFrame:
    rng = np.random.default_rng(cfg.seed)
    bins = list(range(1, cfg.n_bins + 1))
    conc_strs = [str(c) for c in cfg.concs]
    rows: list[dict] = []
    for k in range(cfg.n_clonotypes):
        kd_k = cfg.true_kd * math.exp(rng.normal(0.0, 0.1))
        clonotype = f"C{k:05d}"
        for i, c in enumerate(cfg.concs):
            target = _hill(c, cfg.baseline, cfg.amplitude, kd_k, cfg.true_n)
            centers = np.array(bins, dtype=float)
            probs = np.exp(-0.5 * ((centers - target) / 0.5) ** 2)
            probs /= probs.sum()
            mean_reads = cfg.reads_per_conc * probs
            bin_reads = rng.poisson(mean_reads)
            for j, b in enumerate(bins):
                rows.append(
                    {
                        "clonotypeKey": clonotype,
                        "sampleId": f"s_c{i}_b{b}",
                        "concentrationStr": conc_strs[i],
                        "concentration": c,
                        "bin": b,
                        "reads": int(bin_reads[j]),
                    }
                )
    return pl.DataFrame(rows)


def build_poisson_reads_with_skewed_sort(
    cfg: PoissonConfig,
    sort_fractions_by_conc: dict[float, list[float]],
) -> pl.DataFrame:
    """Skewed-sort variant: cells are preferentially sorted into particular bins.

    Reads at (bin, concentration) are scaled by `sort_fractions_by_conc[conc][bin-1]`
    after Poisson sampling against the Hill-weighted target — this mimics a real
    Tite-Seq sort where the wet-lab gate widths leave some bins much more sparsely
    populated than others. The pipeline's uncorrected formula treats every bin
    equally and so biases mean_bin; the FACS correction divides this bias out.
    """
    rng = np.random.default_rng(cfg.seed)
    bins = list(range(1, cfg.n_bins + 1))
    conc_strs = [str(c) for c in cfg.concs]
    rows: list[dict] = []
    for k in range(cfg.n_clonotypes):
        kd_k = cfg.true_kd * math.exp(rng.normal(0.0, 0.1))
        clonotype = f"C{k:05d}"
        for i, c in enumerate(cfg.concs):
            target = _hill(c, cfg.baseline, cfg.amplitude, kd_k, cfg.true_n)
            centers = np.array(bins, dtype=float)
            probs = np.exp(-0.5 * ((centers - target) / 0.5) ** 2)
            probs /= probs.sum()
            # Scale expected reads per bin by the sort yield at that bin.
            fractions = np.array(sort_fractions_by_conc[c])
            mean_reads = cfg.reads_per_conc * probs * (fractions * cfg.n_bins)
            bin_reads = rng.poisson(mean_reads)
            for j, b in enumerate(bins):
                rows.append(
                    {
                        "clonotypeKey": clonotype,
                        "sampleId": f"s_c{i}_b{b}",
                        "concentrationStr": conc_strs[i],
                        "concentration": c,
                        "bin": b,
                        "reads": int(bin_reads[j]),
                        "sort_fraction": float(fractions[j]),
                    }
                )
    return pl.DataFrame(rows)


@pytest.mark.slow
def test_synthetic_titeseq_facs_correction_improves_recovery():
    """Synthetic skew: uncorrected K_D recovery degrades; corrected recovers it.

    The skew pattern mirrors a real Tite-Seq sort — low-concentration samples
    pile cells into the bottom bin (weak binders stay unlabeled), high-
    concentration samples pile into the top bin (strong binders saturate).
    The uncorrected formula weights every bin equally, so it sees phantom
    signal in the over-represented bins; the corrected formula reweights by
    (C_bc/C_c) and recovers the true Hill shape.

    This test is a deep regression alarm rather than a calibration bound — it
    asserts that the corrected path outperforms the uncorrected path on the
    same fixture, not that either hits a specific numeric threshold.
    """
    cfg = PoissonConfig(
        true_kd=1e-8,
        true_n=1.0,
        baseline=1.0,
        amplitude=math.log(2.5),
        concs=[1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7],
        n_clonotypes=40,
        n_bins=4,
        reads_per_conc=2000,
        seed=7,
    )
    # Heavily skew toward bin 1 at low [antigen] and bin 4 at high [antigen] —
    # the exact failure mode the correction targets.
    sort_fractions_by_conc: dict[float, list[float]] = {
        1e-10: [0.55, 0.25, 0.15, 0.05],
        3e-10: [0.45, 0.30, 0.15, 0.10],
        1e-9: [0.40, 0.30, 0.20, 0.10],
        3e-9: [0.30, 0.30, 0.25, 0.15],
        1e-8: [0.25, 0.25, 0.25, 0.25],
        3e-8: [0.15, 0.25, 0.30, 0.30],
        1e-7: [0.10, 0.20, 0.30, 0.40],
        3e-7: [0.05, 0.15, 0.25, 0.55],
    }
    reads = build_poisson_reads_with_skewed_sort(cfg, sort_fractions_by_conc)

    def _within_factor(out: dict[str, pl.DataFrame], factor: float) -> float:
        pc = out["per_clonotype"]
        good = pc.filter(pl.col("affinityClass").is_in(["Good", "Partial"]))
        if good.height == 0:
            return 0.0
        within = good.filter(
            (pl.col("kd") >= cfg.true_kd / factor) & (pl.col("kd") <= cfg.true_kd * factor)
        )
        return within.height / good.height

    out_corrected = run(reads, sort_fraction_col="sort_fraction")
    out_uncorrected = run(reads)

    corrected_rate = _within_factor(out_corrected, factor=5.0)
    uncorrected_rate = _within_factor(out_uncorrected, factor=5.0)

    # Corrected must recover at least 70% of fits within 5x of the injected K_D.
    # Uncorrected on this skew should be meaningfully worse — assert the gap to
    # catch regressions where the correction silently becomes a no-op.
    assert corrected_rate >= 0.7, (
        f"corrected recovery {corrected_rate:.0%} — correction broken or fixture too hard"
    )
    assert corrected_rate > uncorrected_rate, (
        f"corrected {corrected_rate:.0%} did not beat uncorrected {uncorrected_rate:.0%} — "
        "is the correction actually reaching normalize()?"
    )


@pytest.mark.slow
def test_poisson_titeseq_kd_recovery():
    cfg = PoissonConfig(
        # Sub-µM grid spanning 0.1 nM → 300 nM; K_D placed mid-grid at 10 nM.
        true_kd=1e-8,
        true_n=1.0,
        baseline=1.0,
        amplitude=math.log(2.5),
        concs=[1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7],
        n_clonotypes=50,
        n_bins=4,
        reads_per_conc=2000,
        seed=42,
    )
    reads = build_poisson_reads(cfg)
    out = run(reads)
    pc = out["per_clonotype"]
    good = pc.filter(pl.col("affinityClass").is_in(["Good", "Partial"]))
    # Fraction recovered (not Failed) — deep regression alarm, not a calibration bound.
    frac = good.height / pc.height
    assert frac >= 0.5, f"only {frac:.0%} non-failed (regression alarm, expected >= 50%)"

    # Of those, 50% should have K_D within a factor of 5
    within = good.filter((pl.col("kd") >= cfg.true_kd / 5) & (pl.col("kd") <= cfg.true_kd * 5))
    frac_within = within.height / max(good.height, 1)
    assert frac_within >= 0.5, f"only {frac_within:.0%} within 5x (expected >= 50%)"
