"""Default parameters and shared constants for the titeseq-analysis fitting script."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

# "Below range" and "Above range" are written by "log_fluorescence" only: they say which wall a
# bound K_D sits at. "hill_gate_index" keeps the three classes it shipped with.
AffinityClass = Literal["Good", "Partial", "Below range", "Above range", "Failed"]

# Where the reported K_D sits against the titrated range. Null when the fit set no K_D.
KdRangeSide = Literal["Below range", "In range", "Above range"]

FailureReason = Literal[
    "insufficient_reads",
    "insufficient_points",
    "non_monotonic_signal",
    "convergence_failure",
    "low_r2",
    "high_kd_se",
    "incomplete_titration",
    "n_out_of_range",
    "all_replicates_failed",
    "all_variants_failed",
    "conflicting_censoring",
    "stop_codon",
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
# Starr 2020/2022. They are not in BlockData and not in the UI.
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
# Highest log10 gate fluorescence accepted. Phillips 2021's gate positions run 2.1 to 4.4; linear
# fluorescence routinely exceeds 7. A larger value is linear fluorescence where log10 is required.
MAX_LOG10_FLUORESCENCE: float = 7.0


@dataclass(frozen=True)
class FitParams:
    # The two exposed knobs.
    #
    # This default is deliberately the LEGACY equation, and it is not the default a user gets. The
    # UI sets "log_fluorescence" on a project when the project is created and sends it explicitly.
    # The model omits the key whenever the equation is the legacy one, so an upgraded project
    # projects byte-identical workflow args and does not re-run. An absent key on the wire therefore
    # means "a pre-existing project running v2.4.0's fit", and it must keep meaning that.
    equation: Equation = "hill_gate_index"
    # Read by "log_fluorescence" only. A params file for "hill_gate_index" may omit it.
    plateau_floor: float = 0.667
    # Read by "log_fluorescence" only. An in-range fit with R² below this rates Failed, low_r2.
    # None turns the gate off.
    r2_gate_min: float | None = None
    # Read by "log_fluorescence" only. An in-range fit with SE(log10 K_D) above this rates Failed,
    # high_kd_se. None sets no limit.
    kd_se_gate_max: float | None = None
    # True rates a clonotype with no reads at one or more of the titration's concentrations, 0 M
    # included, Failed, incomplete_titration, before any fit. Phillips 2021 fits complete
    # titrations only.
    require_all_concentrations: bool = False
    # Read by both equations.
    min_reads_per_concentration: int = 3
    min_concentration_points: int = 5
    r2_threshold_good: float = 0.8
    r2_threshold_failed: float = 0.5
    n_min: float = 0.5
    n_max: float = 2.0
    hook_effect_threshold_bin: float = 0.2
    hook_effect_threshold_no_bin: float = 0.02
    hook_effect_min_reads: int = 20

    def __post_init__(self) -> None:
        # `equation` selects the code path, so an unknown value would otherwise run the legacy fit
        # silently. Every other field is range-checked by the stage that reads it.
        if self.equation not in EQUATIONS:
            raise ValueError(f"unknown equation {self.equation!r}; expected one of {', '.join(EQUATIONS)}")
        if self.r2_gate_min is not None and not 0.0 <= self.r2_gate_min <= 1.0:
            raise ValueError(f"r2_gate_min must be between 0 and 1, got {self.r2_gate_min!r}")
        if self.kd_se_gate_max is not None and not self.kd_se_gate_max > 0.0:
            raise ValueError(f"kd_se_gate_max must be above 0, got {self.kd_se_gate_max!r}")


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
# The concentration as a category label that sorts as text in numeric order: "01 · 3.16e-10".
COL_CONC_LABEL = "concentrationLabel"
COL_BIN = "bin"
COL_ANTIGEN = "antigen"
COL_READS = "reads"
# Column name emitted by the Tengo workflow in reads.tsv when
# sortFractionColumnRef is bound. Must match the string passed to
# --sort-fraction-column. Do not change without a simultaneous change
# to workflow/src/main.tpl.tengo.
COL_SORT_FRACTION = "sort_fraction"
# Column header the Tengo workflow writes into reads.tsv when binFluorescenceColumnRef is bound:
# the gate's measured mean log10 fluorescence, one value per (concentration, gate). The pipeline
# detects the column by this name; no CLI flag carries it. Do not change without a simultaneous
# change to workflow/src/main.tpl.tengo. Read by "log_fluorescence" only. Absent means the gates
# are assumed evenly spaced at GATE_WIDTH_DECADES.
COL_BIN_FLUORESCENCE = "bin_fluorescence"
# Column header the Tengo workflow writes into reads.tsv when replicateColumnRef is bound: which
# sample replicate each sample belongs to. Named on the command line by --replicate-column. Do not
# change without a simultaneous change to workflow/src/main.tpl.tengo.
COL_REPLICATE = "replicate"
# The label every row carries when no replicate column is bound. The replicate combination does not
# run then, and no output carries a replicate.
SINGLE_REPLICATE = "1"
