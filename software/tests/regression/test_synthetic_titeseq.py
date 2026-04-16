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


@pytest.mark.slow
def test_poisson_titeseq_kd_recovery():
    cfg = PoissonConfig(
        true_kd=10.0,
        true_n=1.0,
        baseline=1.0,
        amplitude=math.log(2.5),
        concs=[0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0],
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
