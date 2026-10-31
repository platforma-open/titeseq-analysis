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

# The equation fitted to each clonotype.
#   "log_fluorescence"  — log10(A*c/(c+K) + B) on mean log10 gate fluorescence, exponent fixed at 1,
#                         weighted least squares. Requires gates (bin mode).
#   "hill_gate_index"   — baseline + (top-baseline)*c^n/(K^n+c^n) on the read-weighted mean gate
#                         index, exponent free. The v2.4.0 equation; also runs in no-bin mode.
Equation = Literal["log_fluorescence", "hill_gate_index"]
EQUATIONS: tuple[Equation, ...] = ("log_fluorescence", "hill_gate_index")

# Minimum dynamic range between top and baseline (R10).
# Bin mode: δ = 0.5 on bin-index scale; no-bin mode: δ = 0.05 on frequency scale.
DELTA_BIN: float = 0.5
DELTA_NO_BIN: float = 0.05

# --- Fixed internal constants for "log_fluorescence". ------------------------------------------
# These are constants, not parameters: each takes one value across Adams 2016, Phillips 2021 and
# Starr 2020/2022. They are not in BlockData and not in the UI. The equation ignores the user's
# min_reads_per_concentration / min_concentration_points and uses the two below instead.
LF_MIN_READS_PER_CONCENTRATION: int = 3
LF_MIN_CONCENTRATION_POINTS: int = 5
# Observed rise, max(F̄) − min(F̄), at or above which a curve counts as a measurement rather than a
# bound. Feeds the rating's bound test only; it changes no fitted K_D.
LF_RISE_THRESHOLD_DECADES: float = 0.70
# Standard error of log10 K_D at or below which an in-range fit is rated Good.
KD_SD_GOOD: float = 0.5
# Width of one FACS gate in decades of fluorescence. Gate label g sits at (g − 1) * this.
GATE_WIDTH_DECADES: float = 0.65
# Gate range assumed when the gate labels do not define one (a single gate label).
GATE_RANGE_DECADES: float = 1.950
# The exponent "log_fluorescence" fixes.
HILL_EXPONENT_FIXED: float = 1.0


@dataclass(frozen=True)
class FitParams:
    # The two exposed knobs.
    #
    # This default is deliberately the LEGACY equation, and it is not the default a user gets. The
    # UI defaults a new project to "log_fluorescence" and sends it explicitly. The model omits the
    # key whenever the equation is the legacy one, so an upgraded project projects byte-identical
    # workflow args and does not re-run. An absent key on the wire therefore means "a pre-existing
    # project running v2.4.0's fit", and it must keep meaning that.
    equation: Equation = "hill_gate_index"
    # Read by "log_fluorescence" only. A params file for "hill_gate_index" may omit it.
    plateau_floor: float = 2.0 / 3.0
    # Everything below is read by "hill_gate_index" only.
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
# Canonical concentration string preserved byte-for-byte through the pipeline
# (R14). Acts as the join axis key in output PColumns. The Tengo workflow wraps
# this column as a `pl7.app/vdj/concentration` axis of type String.
COL_CONC_STR = "concentrationStr"
# Internal Float64 used for arithmetic (Hill fit, baseline, weights). Never
# written into output TSVs as the axis key — output rows carry COL_CONC_STR
# only, so no float→string→float roundtrip happens between Python and the
# Tengo xsv.importFile boundary.
COL_CONC_VAL = "concentration"
COL_BIN = "bin"
COL_ANTIGEN = "antigen"
COL_READS = "reads"
# Column name emitted by the Tengo workflow in reads.tsv when
# sortFractionColumnRef is bound. Must match the string passed to
# --sort-fraction-column. Do not change without a simultaneous change
# to workflow/src/main.tpl.tengo.
COL_SORT_FRACTION = "sort_fraction"
