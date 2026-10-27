"""Synthetic-data builders shared across the test suite. Independent of src: these are the reference."""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from constants import FitParams

# The opt-in linear Hill model on mean bin. Fixtures without cell counts or a 0 M baseline need it.
LINEAR_HILL = FitParams(fit_model="hill_free_n")


def hill_curve(x: np.ndarray, baseline: float, amplitude: float, kd: float, n: float) -> np.ndarray:
    """Hill curve with top = baseline + e^amplitude. A reference for the fitter, so it does not reuse src."""
    top = baseline + math.exp(amplitude)
    return baseline + (top - baseline) * (x**n) / (kd**n + x**n)


OUTPUTS = ("per_clonotype", "mean_bin", "fitted_mean_bin", "concentration_value")


def cli_outputs(out_dir: Path, suffix: str = ".tsv") -> tuple[list[str], dict[str, Path]]:
    """The four required output arguments of main.py, writing into out_dir, and the paths they name."""
    paths = {name: out_dir / f"{name}{suffix}" for name in OUTPUTS}
    argv = [x for name, path in paths.items() for x in (f"--out-{name.replace('_', '-')}", str(path))]
    return argv, paths


def hill_reads(
    clonotype: str,
    true_kd: float,
    concs: list[float],
    bins: list[int],
    *,
    baseline: float = 1.5,
    amplitude: float = math.log(2.0),
    n: float = 1.0,
    per_conc: int = 500,
    cells: Callable[[int], float] | None = None,
    antigen: str | None = None,
) -> list[dict]:
    """Noiseless reads for one clonotype: per_conc reads per concentration, spread over the bins as a Gaussian
    (sd 0.35 bins) around the Hill mean bin `hill_curve(c, baseline, amplitude, true_kd, n)`.

    cells(bin) adds a `cells` column; antigen adds an `antigen` column and prefixes the sample ids.
    """
    targets = hill_curve(np.asarray(concs, dtype=float), baseline, amplitude, true_kd, n)
    prefix = f"{antigen}_" if antigen is not None else ""
    rows = []
    for i, c in enumerate(concs):
        weights = np.exp(-0.5 * ((np.array(bins, dtype=float) - targets[i]) / 0.35) ** 2)
        weights /= weights.sum()
        counts = np.round(weights * per_conc).astype(int)
        if counts.sum() == 0:
            counts[len(counts) // 2] = 1
        for j, b in enumerate(bins):
            row = {
                "clonotypeKey": clonotype,
                "sampleId": f"{prefix}s_c{i}_b{b}",
                "concentrationStr": str(c),
                "concentration": float(c),
                "bin": int(b),
                "reads": int(counts[j]),
            }
            if cells is not None:
                row["cells"] = cells(b)
            if antigen is not None:
                row["antigen"] = antigen
            rows.append(row)
    return rows


@dataclass(frozen=True)
class SimulatedExperiment:
    reads: pl.DataFrame  # block input, with a `cells` column (cells sorted into the sample)
    true_log10_kd: dict[str, float]
    bin_width: float
    baseline: float


def simulate_titeseq(
    clonotypes: dict[str, tuple[float, float]],
    *,
    concentrations: Sequence[float] = (0.0,) + tuple(10 ** (-10 + 0.5 * i) for i in range(10)),
    n_bins: int = 4,
    bin_width: float = 0.7,
    baseline: float = 2.5,
    cell_sd: float = 0.3,
    cells_per_clonotype: int = 4000,
    reads_per_sample: float = 2e5,
    seed: int = 0,
) -> SimulatedExperiment:
    """Tite-Seq at read level, after Adams et al. 2016 Appendix G.

    clonotypes: key -> (log10 K_D, amplitude A in bin-1 fluorescence units). Each cell's log10
    fluorescence is Normal(log10(B + A c/(c+K)), cell_sd); FACS puts it in the bin whose value
    bin_width * b is nearest. Sequencing draws Poisson reads in proportion to each clonotype's share
    of the bin's cells.
    """
    rng = np.random.default_rng(seed)
    keys = sorted(clonotypes)
    rows = []
    for ci, c in enumerate(concentrations):
        counts = np.zeros((len(keys), n_bins))
        for i, key in enumerate(keys):
            log_kd, amp = clonotypes[key]
            mu = np.log10(baseline + amp * c / (c + 10.0**log_kd))
            x = rng.normal(mu, cell_sd, cells_per_clonotype)
            b = np.clip(np.rint(x / bin_width), 0, n_bins - 1).astype(int)
            counts[i] = np.bincount(b, minlength=n_bins)
        cells = counts.sum(0)
        for b in range(n_bins):
            if cells[b] == 0:
                continue
            reads = rng.poisson(reads_per_sample * counts[:, b] / cells[b])
            for i, key in enumerate(keys):
                if reads[i] > 0:
                    rows.append(
                        {
                            "clonotypeKey": key,
                            "sampleId": f"c{ci}_b{b + 1}",
                            "concentrationStr": repr(c) if c else "0",
                            "concentration": float(c),
                            "bin": b + 1,
                            "reads": int(reads[i]),
                            "cells": float(cells[b]),
                        }
                    )
    return SimulatedExperiment(pl.DataFrame(rows), {k: v[0] for k, v in clonotypes.items()}, bin_width, baseline)
