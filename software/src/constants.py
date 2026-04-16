"""Default parameters and shared constants for the titeseq-analysis fitting script."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

AffinityClass = Literal["Good", "Partial", "Failed"]

FailureReason = Literal[
    "insufficient_reads",
    "insufficient_points",
    "non_monotonic_signal",
    "convergence_failure",
    "low_r2",
    "n_out_of_range",
]

# Minimum dynamic range between top and baseline (R10).
DELTA: float = 0.05


@dataclass(frozen=True)
class FitParams:
    min_reads_per_concentration: int = 3
    min_concentration_points: int = 5
    r2_threshold_good: float = 0.8
    r2_threshold_failed: float = 0.5
    n_min: float = 0.5
    n_max: float = 2.0
    hook_effect_threshold_bin: float = 0.2
    hook_effect_threshold_no_bin: float = 0.02
    hook_effect_min_reads: int = 20


DEFAULT_PARAMS = FitParams()

# Column name conventions used inside the pipeline.
COL_CLONOTYPE = "clonotypeKey"
COL_SAMPLE = "sampleId"
COL_CONC_STR = "concentrationStr"  # canonical string key (R14)
COL_CONC_VAL = "concentration"  # numeric value
COL_BIN = "bin"
COL_ANTIGEN = "antigen"
COL_READS = "reads"
