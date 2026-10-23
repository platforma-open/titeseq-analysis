"""Shared fixtures."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
import polars as pl
import pytest


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


@pytest.fixture(scope="session")
def simulate() -> Callable[..., SimulatedExperiment]:
    """The read-level simulator, for tests that need their own clonotype set."""
    return simulate_titeseq


@pytest.fixture(scope="session")
def simulated_library() -> SimulatedExperiment:
    """Binders across the titration range, one far tighter, one non-binder, at mixed expression."""
    clonotypes = {}
    for j, log_kd in enumerate([-9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0]):
        for r, amp in enumerate([40.0, 80.0]):
            clonotypes[f"B{j}{r}"] = (log_kd, amp)
    clonotypes["TIGHT"] = (-12.0, 60.0)
    clonotypes["NONBINDER"] = (-2.0, 60.0)
    return simulate_titeseq(clonotypes, seed=1)
