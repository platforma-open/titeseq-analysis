"""Shared fixtures."""

from __future__ import annotations

import pytest
from helpers import SimulatedExperiment, simulate_titeseq


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
