"""Phillips, Lawrence, Moulana et al. (2021), "Binding affinity landscapes constrain the evolution of broadly
neutralizing anti-influenza antibodies", eLife 10:e71393 (CC-BY 4.0). Antibody CR6261.

Values below are copied from the paper text and the authors' scripts. The data files (the authors' read
counts, gate positions, cell counts, their per-variant affinities and the isogenic titrations) are read
from a local copy only, never committed.

The library holds all 2^11 combinations of 11 heavy-chain mutations. A variant is an 11-bit string, bit i
= 1 when mutation i is in its mature state: `00000000000` is the germline, `11111111111` the mature
antibody's 11 library mutations.
"""

from __future__ import annotations

from dataclasses import dataclass

STUDY = "phillips2021-cr6261"

# Two antigens, three biological replicates each; each library is one campaign's titration.
LIBRARIES = {"h1": ("a", "b", "c"), "h9": ("x", "y", "z")}

# The concentration code in the authors' column names (`correct_counts_6261_<lib>_<code>_<bin>`), and in
# the gate and cell-count tables (`_<code>_`), is -log10 of the molar concentration without its point:
# "75" is 10^-7.5 M. "0" is the no-antigen control. Checked against upload/metadata.csv, whose sample
# `cr6261-h1-a-c75-bin1` carries 3.162e-08 M.
NEG_LOG10_CONC: dict[str, float | None] = {
    "0": None,
    "7": 7.0,
    "75": 7.5,
    "8": 8.0,
    "85": 8.5,
    "9": 9.0,
    "95": 9.5,
    "10": 10.0,
    "105": 10.5,
    "11": 11.0,
    "115": 11.5,
    "12": 12.0,
}
N_BINS = 4

# Methods: "half-log increments spanning 1 pM – 100 nM". The authors' fitting script (`Mean bin Kd
# inference with error-6261a.ipynb`) uses 10^-7.05 ... 10^-12.05 M instead, with no reason given, so
# their K_D values read 0.05 decades tighter than the labels give. The paper case reproduces that.
AUTHOR_CONC_SHIFT = 0.05

# Methods: K_D values outside the titration are pinned to 10^-7 and 10^-12 M. Published values are
# -log10 K_D; the unfiltered tables used here are not pinned (per-replicate values run 5.0 to 10.6).
NEG_LOG10_KD_LOWER = 7.0
NEG_LOG10_KD_UPPER = 12.0

GERMLINE = "00000000000"
MATURE = "11111111111"


@dataclass(frozen=True)
class Isogenic:
    """A variant the authors also built as a single clone and titrated by flow cytometry
    (answer-key/isogenic/Isogenic_Kd/Mutants/fitted_Kds.csv, two replicate titrations each).
    `name` is the file's `shortname`: the mutations, numbered from 1, the clone carries."""

    name: str
    genotype: str


# Every clone with a fitted isogenic K_D. "som" is the mature antibody (the 11 library mutations, the
# library's all-ones corner). The Revertants file's "som14" carries all 14 somatic mutations, three of
# them outside the library, so it has no library genotype and is left out. The non-binding clones in
# mut-list.csv (1, 2, 6, germline) have no fitted isogenic K_D.
ISOGENIC = [
    Isogenic("1,2", "11000000000"),
    Isogenic("6,7", "00000110000"),
    Isogenic("1,2,6,7", "11000110000"),
    Isogenic("1,2,6,7,8", "11000111000"),
    Isogenic("som", MATURE),
]
