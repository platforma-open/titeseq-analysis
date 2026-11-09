"""Adams, Mora, Walczak, Kinney (2016), "Measuring the sequence-affinity landscape of antibodies with
massively parallel titration curves", eLife 5:e23156 (CC-BY 4.0).

Values below are copied from the paper text. The authors' data files carry no licence: they are read from
a local copy or downloaded at run time from the repository at PINNED_COMMIT, never committed.
"""

from __future__ import annotations

from dataclasses import dataclass

STUDY = "adams2016"
REPO = "jbkinney/16_titeseq"
PINNED_COMMIT = "54130c8fc6081bbfa8301c41e1afcfa6a5d79590"
RAW_URL = f"https://raw.githubusercontent.com/{REPO}/{PINNED_COMMIT}/data/{{name}}"

# Results B: eleven fluorescein concentrations, 4 sorting bins, 3 replicates. Column labels in
# data/replicate_N.csv are `fluorescein<label>bin<b>`, b from 0.
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
SORT_COUNT_FILES = {1: "sort_counts_16.4.15.txt", 2: "sort_counts_16.4.19.txt", 3: "sort_counts_16.4.21.txt"}

# Fig 4 caption and Table I: K_D is measurable between these limits; values beyond them are bounds.
# The authors' own combine step clips each replicate to them (figure_4_accuracy.py:21-22, 288-310).
LOG10_KD_LOWER = -9.5
LOG10_KD_UPPER = -5.0


@dataclass(frozen=True)
class Clone:
    """Table I row: one nucleotide sequence, measured by Tite-Seq (mean and SE over the 3 replicates)
    and by flow cytometry. `censor` marks a Tite-Seq value at a limit."""

    name: str  # `name` in data/titration_curves.csv, which gives the nucleotide CDRs
    flow_log10_kd: float
    tite_seq_log10_kd: float
    tite_seq_se: float | None
    censor: str | None = None


# Table I, ordered by flow K_D. The published Tite-Seq column cannot be rebuilt exactly from the
# deposited replicate files (the authors' own code gives slightly different numbers), so the test reads
# it as a reference, not as an identity.
TABLE_I = [
    Clone("opt", -9.5, -9.5, None, "below"),
    Clone("107", -9.28, -9.18, 0.11),
    Clone("112", -8.95, -9.19, 0.14),
    Clone("pJK36", -8.61, -8.92, 0.10),
    Clone("144", -8.57, -8.86, 0.04),
    Clone("133", -8.55, -8.62, 0.09),
    Clone("132", -8.48, -8.38, 0.29),
    Clone("94", -8.46, -8.50, 0.04),
    Clone("5", -8.34, -8.55, 0.09),
    Clone("93", -7.35, -7.60, 0.70),
    Clone("39", -7.08, -7.28, 0.17),
    Clone("102", -5.76, -7.25, 0.60),
    Clone("22", -5.69, -7.53, 0.07),
    Clone("7", -5.53, -5.39, 0.18),
    Clone("45", -5.40, -5.0, None, "above"),
    Clone("103", -5.15, -5.44, 0.55),
    Clone("3", -5.0, -5.0, None, "above"),
    Clone("18", -5.0, -5.0, None, "above"),
]
