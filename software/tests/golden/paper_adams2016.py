"""Adams, Mora, Walczak, Kinney (2016), "Measuring the sequence-affinity landscape of antibodies with
massively parallel titration curves", eLife 5:e23156 (CC-BY 4.0).

Values below are copied from the paper text. The authors' data files carry no licence: they are read from
a local copy or downloaded at run time from the repository at PINNED_COMMIT, never committed.
"""

from __future__ import annotations

import math
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

# The FACS recordings of each sort, one folder a replicate: `data/fcs<rep>/<prefix>_Sort <n>_<file no>.fcs`.
# Sort 1 is the expression sort, Sort 2 the 0 M sort, Sorts 3 to 12 the concentrations in ascending
# order (figure_3_procedure.py:183, 220-277).
FCS_PREFIX = {1: "040616", 2: "041916", 3: "04,2f,21,2f,16"}


def fcs_name(rep: int, concentration_index: int) -> str:
    """The recording of one affinity sort, by the index of its concentration in CONCENTRATION_LABELS."""
    sort = concentration_index + 2
    return f"fcs{rep}/{FCS_PREFIX[rep]}_Sort {sort}_{sort + 1:03d}.fcs"


# The authors' gating, in log10 units (figure_3_procedure.py). Cells pass the FSC-A/SSC-A polygon
# (lines 92-98), then the FITC-A/PE-A pre-sort polygon (lines 108-116); the four affinity gates are
# PE-A intervals, set anew on each sort day (lines 184-199).
FSC_SSC_POLYGON = [
    (4.33487278721722, 4.336993061421767),
    (4.004403501790057, 4.010254050207731),
    (3.5539963164652484, 3.324117555127123),
    (3.761167625032098, 3.10094159611305),
    (4.797024167866345, 4.285490917033903),
    (4.701406640835492, 4.474332113122737),
    (4.33487278721722, 4.336993061421767),
]
PRE_SORT_POLYGON = [
    (2.1356696655075864, 1.6073794088650155),
    (2.8528011182389887, 1.6588815532528791),
    (2.9802911542801267, 2.6717570595475224),
    (4.637661622814923, 4.491499494585357),
    (3.7133588615166713, 5.058023082851852),
    (1.7850720663944564, 5.040855701389231),
    (1.5141557398070378, 1.7790532234912266),
    (1.6575820303533182, 1.487207738626668),
    (2.1356696655075864, 1.6073794088650155),
]
LOG10_PE_GATES = {
    1: [
        (1.4775788014225415, 2.245106453825187),
        (2.2459637507318497, 2.9683536647430646),
        (2.96837636283699, 3.6663288893768895),
        (3.6675170811229103, math.log10(3e4)),
    ],
    2: [
        (1.4775788014225415, 2.245106453825187),
        (2.2459637507318497, 2.85027248745035),
        (2.8650326346119392, 3.474446976276228),
        (3.4849588153986915, math.log10(1e5)),
    ],
    3: [
        (1.4775788014225415, 2.200826012340419),
        (2.216443456408671, 2.8355123402887603),
        (2.8355123402887603, 3.474446976276228),
        (3.4849588153986915, math.log10(3e4)),
    ],
}

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
