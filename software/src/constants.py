"""Default parameters and shared constants for the titeseq-analysis fitting script."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

# "CensoredBelow": K_D sits at or under the lowest titrated concentration, so the reported value is
# an upper bound on affinity, not a measurement. Its own class, outside the ["Good","Partial"]
# default cutoff, because the bound value is the smallest in the table and would otherwise sort to
# the top of any downstream ranking as an unorderable block.
AffinityClass = Literal["Good", "Partial", "CensoredBelow", "Failed"]

# A converged fit is never Failed. Every reason below is a refusal to fit at all.
FailureReason = Literal[
    "insufficient_reads",
    "insufficient_points",
    "non_monotonic_signal",
    "convergence_failure",
    "insufficient_dynamic_range",
]

# Minimum dynamic range between top and baseline (R10).
# Bin mode: δ = 0.5 on bin-index scale; no-bin mode: δ = 0.05 on frequency scale.
DELTA_BIN: float = 0.5
DELTA_NO_BIN: float = 0.05


# Curve model for bin mode.
#   "poisson_log": Adams, Mora, Walczak, Kinney 2016 Appendix E — Poisson likelihood on read counts,
#                  mean log fluorescence = log10(B + A*c/(c+K)), Hill n = 1. The default.
#   "hill_free_n": linear Hill with free n fitted to mean bin by weighted least squares. Opt-in;
#                  does not reproduce Adams 2016.
FitModel = Literal["poisson_log", "hill_free_n"]


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
    # R9b hook detector. Off by default: with the Hill exponent pinned at 1 the fits it flags fail to
    # converge instead, and leaving it on costs 44/99/73 of the 250 tightest binders on the three
    # Adams 2016 replicates while buying no rank agreement.
    hook_detector_enabled: bool = False
    # Fix the Hill exponent at 1, as the Adams 2016 model does. Free n costs rank agreement
    # (0.909 -> 0.816 on Adams replicate 1) and is the condition under which the hook detector is
    # redundant. Turn off to fit n.
    pin_hill_exponent: bool = True
    fit_model: FitModel = "poisson_log"
    # Width of one FACS bin in decades of fluorescence (gate spacing). Adams 2016: 0.61-0.71.
    bin_width_decades: float = 0.65
    # Amplitude prior: every clonotype must be able to reach this fraction of the bin range at
    # saturation. Adams 2016 used 0.63-0.70. Flat curves then fit as K_D above the titration range.
    min_plateau_fraction: float = 0.65
    # Mean fluorescence at 0 M in units of bin-1 fluorescence. None: computed from the 0 M samples.
    baseline_fluorescence: float | None = None
    # Amplitude floor as a fraction of the bin range, applied inside the Hill fit bounds. None keeps
    # the floor at DELTA_BIN. Raising it deletes rows whose plateau is low for reasons unrelated to
    # affinity, so it trades retention for the average quality of what remains.
    amplitude_prior_q: float | None = None
    # Posterior SD of log10 K_D at or below which an in-range K_D is classed Good.
    kd_sd_good: float = 0.5


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
# Column name for cells sorted into the sample (absolute count). Same contract as COL_SORT_FRACTION:
# must match the string passed to --cell-count-column.
COL_CELL_COUNT = "cell_count"
# Column name for the mean log10 fluorescence of the cells sorted into the sample. Same contract as
# COL_SORT_FRACTION: must match the string passed to --bin-fluorescence-column.
COL_BIN_FLUORESCENCE = "bin_fluorescence"
