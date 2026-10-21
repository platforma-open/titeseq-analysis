"""Golden values from Adams, Mora, Walczak, Kinney (2016) "Measuring the sequence-affinity landscape of
antibodies with massively parallel titration curves", eLife 5:e23156 (CC-BY 4.0).

Numbers here are copied from the paper text. Where a value is an analysis setting the paper does not
print, it is derived from the authors' code at PINNED_COMMIT and the derivation is stated next to it.
The authors' data files carry no licence: tests download them at run time and never commit them.
"""

from __future__ import annotations

from dataclasses import dataclass

REPO = "jbkinney/16_titeseq"
PINNED_COMMIT = "54130c8fc6081bbfa8301c41e1afcfa6a5d79590"

# Main text, Results B: eleven fluorescein concentrations, 4 sorting bins, 3 replicates.
CONCENTRATION_LABELS = [
    "0",
    "10^-9.5",
    "10^-9",
    "10^-8.5",
    "10^-8",
    "10^-7.5",
    "10^-7",
    "10^-6.5",
    "10^-6",
    "10^-5.5",
    "10^-5",
]
CONCENTRATIONS_M = [0.0] + [10 ** (-9.5 + 0.5 * i) for i in range(10)]
N_BINS = 4

# Fig 4 caption and Table I: K_D is measurable between these limits; values beyond them are bounds.
LOG10_KD_LOWER_LIMIT = -9.5
LOG10_KD_UPPER_LIMIT = -5.0

# Results D: WT averaged over 33 synonymous variants and all three replicates.
WT_CDR1H_AA = "TFSDYWMNWV"
WT_CDR3H_AA = "GSYYGMDYWG"
WT_LOG10_KD = -8.87
WT_LOG10_KD_SE = 0.02

# Fig 4C: coefficient of determination between log Tite-Seq and log flow K_D over the Table I clones.
FIG4C_R2 = 0.89


@dataclass(frozen=True)
class Replicate:
    number: int
    sort_counts_file: str
    # Gate width in decades of fluorescence. The paper does not print it. Derived two ways from the
    # authors' code, which agree: the bin gates in validation_simulations.py (bin_vals_16_4_*), and
    # log10(max fit_saturation)/3, the top-bin value of the amplitude grid in replicate_N.csv.
    bin_width_decades: float
    # Lower bound of the amplitude grid, expressed as the minimum saturated plateau as a fraction of
    # the bin range: log10(B + min fit_saturation) / bin_width / 3, with B the replicate's `basal`.
    # The paper does not state how the authors chose it.
    min_plateau_fraction: float
    # B, the 0 M mean fluorescence in bin-1 units: the `basal` column of replicate_N.csv. The block
    # computes it from the 0 M reads (golden test); the fast tier passes it to skip the full library.
    baseline_fluorescence: float


REPLICATES = [
    Replicate(1, "sort_counts_16.4.15.txt", 0.71017, 0.6427, 2.59256067),
    Replicate(2, "sort_counts_16.4.19.txt", 0.60686, 0.6318, 4.23348966),
    Replicate(3, "sort_counts_16.4.21.txt", 0.62900, 0.6950, 2.69899748),
]


@dataclass(frozen=True)
class Clone:
    """Table I row. `censor` is "below" for values at or under 10^-9.5 and "above" for 10^-5.0 or over."""

    name: str  # name in data/titration_curves.csv
    cdr1h_aa: str
    cdr3h_aa: str
    flow_log10_kd: float
    tite_seq_log10_kd: float
    tite_seq_se: float | None
    tite_seq_censor: str | None = None


# Table I, ordered by flow K_D. Tite-Seq values are the mean and SE across the three replicates
# for one specific nucleotide variant (not averaged over synonymous variants).
TABLE_I = [
    Clone("opt", "TFGHYWMNWV", "GASYGMEYLG", -9.5, -9.5, None, "below"),
    Clone("107", "TFSDYWMNWV", "GAYYGMDYWG", -9.28, -9.18, 0.11),
    Clone("112", "TFSDYWMNWV", "GSYYGMDYCG", -8.95, -9.19, 0.14),
    Clone("pJK36", "TFSDYWMNWV", "GSYYGMDYWG", -8.61, -8.92, 0.10),
    Clone("144", "VFSDYWMNWV", "GSYYGMDYWG", -8.57, -8.86, 0.04),
    Clone("133", "AFSDYWMNWV", "GSYYGMDYWG", -8.55, -8.62, 0.09),
    Clone("132", "TFMDYWLNWV", "GSYYGMDYWG", -8.48, -8.38, 0.29),
    Clone("94", "TFSDYWMNWV", "GSYYGMDSWG", -8.46, -8.50, 0.04),
    Clone("5", "TFSDYWINWV", "GSYYGMDYWG", -8.34, -8.55, 0.09),
    Clone("93", "TFSDYWMNWV", "GSYRGMDYWG", -7.35, -7.60, 0.70),
    Clone("39", "TFSDYWMNWV", "GSYYGMDYWA", -7.08, -7.28, 0.17),
    Clone("102", "TFSDYWMNWV", "SSKYGMDYWG", -5.76, -7.25, 0.60),
    Clone("22", "SSSDYWMNWV", "GSYYGMDYWG", -5.69, -7.53, 0.07),
    Clone("7", "HFSDYWMNWL", "GSYYGMDYWG", -5.53, -5.39, 0.18),
    Clone("45", "TFSDYWMNWV", "GSYDGNDYWG", -5.40, -5.0, None, "above"),
    Clone("103", "TFSDYWMNWV", "GSYYGMDLWG", -5.15, -5.44, 0.55),
    Clone("3", "TFSDYWMSWV", "GSYYGMDYWG", -5.0, -5.0, None, "above"),
    Clone("18", "TFSDYSMNWV", "GSYYGMDYWG", -5.0, -5.0, None, "above"),
    Clone("delta", "", "", -5.0, -5.0, None, "above"),
]
