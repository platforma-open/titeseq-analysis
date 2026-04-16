"""Shared fixtures and data builders."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import polars as pl
import pytest

from constants import DEFAULT_PARAMS, FitParams


@pytest.fixture
def default_params() -> FitParams:
    return DEFAULT_PARAMS


def _rows(clonotype: str, conc_strs: list[str], bins: list[int] | None, reads_matrix: list[list[int]],
          antigen: str | None = None, sample_prefix: str | None = None) -> list[dict]:
    """Build long-form rows for a single clonotype.

    If bins is None => no-bin mode; reads_matrix is list[list[int]] with one row per concentration
    but only first element used. If bins provided, reads_matrix[i][j] = reads for conc i, bin j.
    """
    out = []
    for i, cs in enumerate(conc_strs):
        if bins is None:
            row = {
                "clonotypeKey": clonotype,
                "sampleId": f"{sample_prefix or clonotype}_c{i}",
                "concentrationStr": cs,
                "concentration": float(cs),
                "reads": int(reads_matrix[i][0]),
            }
            if antigen is not None:
                row["antigen"] = antigen
            out.append(row)
        else:
            for j, b in enumerate(bins):
                row = {
                    "clonotypeKey": clonotype,
                    "sampleId": f"{sample_prefix or clonotype}_c{i}_b{b}",
                    "concentrationStr": cs,
                    "concentration": float(cs),
                    "bin": int(b),
                    "reads": int(reads_matrix[i][j]),
                }
                if antigen is not None:
                    row["antigen"] = antigen
                out.append(row)
    return out


def build_reads(
    clonotypes: list[str],
    conc_strs: list[str],
    bins: list[int] | None,
    reads_by_clonotype: dict[str, list[list[int]]],
    antigens: dict[str, str] | None = None,
) -> pl.DataFrame:
    """Construct a canonical long-format reads frame for tests.

    `reads_by_clonotype[clonotype]` is list[list[int]]: outer index = concentration index,
    inner = bin index (or [reads_total] when bins is None).
    """
    all_rows: list[dict] = []
    # One sample per concentration needs shared reads across clonotypes, so use
    # a consistent sample scheme: same sampleId for same (concentration, bin) across clonotypes.
    for clonotype in clonotypes:
        antigen = antigens.get(clonotype) if antigens else None
        for i, cs in enumerate(conc_strs):
            if bins is None:
                all_rows.append(
                    {
                        "clonotypeKey": clonotype,
                        "sampleId": f"s_c{i}",
                        "concentrationStr": cs,
                        "concentration": float(cs),
                        "reads": int(reads_by_clonotype[clonotype][i][0]),
                        **({"antigen": antigen} if antigen else {}),
                    }
                )
            else:
                for j, b in enumerate(bins):
                    all_rows.append(
                        {
                            "clonotypeKey": clonotype,
                            "sampleId": f"s_c{i}_b{b}",
                            "concentrationStr": cs,
                            "concentration": float(cs),
                            "bin": int(b),
                            "reads": int(reads_by_clonotype[clonotype][i][j]),
                            **({"antigen": antigen} if antigen else {}),
                        }
                    )
    return pl.DataFrame(all_rows)


@dataclass
class PoissonConfig:
    true_kd: float
    true_n: float
    baseline: float
    amplitude: float  # log(top - baseline)
    concs: list[float]
    n_clonotypes: int
    n_bins: int
    reads_per_conc: int
    seed: int


def _hill_noise_free(x: float, baseline: float, amplitude: float, kd: float, n: float) -> float:
    top = baseline + math.exp(amplitude)
    return baseline + (top - baseline) * (x**n) / (kd**n + x**n)


def build_poisson_reads(cfg: PoissonConfig) -> pl.DataFrame:
    """Simulate reads at the bin level. For each clonotype and concentration, the expected
    mean bin is set by the Hill curve; we draw read counts per bin from a Poisson distribution
    whose per-bin rate reproduces that mean bin (approximately).
    """
    rng = np.random.default_rng(cfg.seed)
    bins = list(range(1, cfg.n_bins + 1))
    conc_strs = [str(c) for c in cfg.concs]
    rows: list[dict] = []
    for k in range(cfg.n_clonotypes):
        # Add per-clonotype jitter on true K_D to widen the regression target spread.
        kd_k = cfg.true_kd * math.exp(rng.normal(0.0, 0.1))
        clonotype = f"C{k:05d}"
        for i, c in enumerate(cfg.concs):
            target_mean = _hill_noise_free(c, cfg.baseline, cfg.amplitude, kd_k, cfg.true_n)
            # Distribute reads so that sum_{b}(b * reads_b / depth_b) / sum = target_mean.
            # Simpler: bias probability of bin b by a Gaussian centered at target_mean.
            centers = np.array(bins, dtype=float)
            probs = np.exp(-0.5 * ((centers - target_mean) / 0.5) ** 2)
            probs /= probs.sum()
            mean_reads_per_bin = cfg.reads_per_conc * probs
            bin_reads = rng.poisson(mean_reads_per_bin)
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
