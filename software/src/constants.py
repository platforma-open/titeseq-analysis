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
# Bin mode: δ = 0.5 on bin-index scale; no-bin mode: δ = 0.05 on frequency scale.
DELTA_BIN: float = 0.5
DELTA_NO_BIN: float = 0.05


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
COL_CONC_STR = "concentrationStr"  # canonical string key (R14) — internal joins only
COL_CONC_VAL = "concentration"  # numeric value (assumed Molar)
COL_CONC_AM = "concentrationAM"  # attomolar Int64 — axis key in output TSVs
CONC_AM_SCALE = 1_000_000_000_000_000_000  # 1e18: Molar → attomolar
# Attomolar encoding uses Int64; anything above this overflows. TiteSeq is a sub-µM
# assay, so this ceiling (~9.22 M) is well above any realistic input. A value above
# the ceiling almost certainly indicates a unit-entry mistake (e.g., molar vs molal,
# or typed "100" instead of "1e-7").
MAX_CONCENTRATION_M = (2**63 - 1) / CONC_AM_SCALE
COL_BIN = "bin"
COL_ANTIGEN = "antigen"
COL_READS = "reads"
# Column name emitted by the Tengo workflow in reads.tsv when
# sortFractionColumnRef is bound. Must match the string passed to
# --sort-fraction-column. Do not change without a simultaneous change
# to workflow/src/main.tpl.tengo.
COL_SORT_FRACTION = "sort_fraction"
