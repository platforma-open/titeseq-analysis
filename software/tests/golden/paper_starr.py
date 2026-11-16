"""Starr 2020 and Starr 2022: what the papers published, and how their analysis code fits.

Both are yeast-display RBD libraries titrated against hACE2: four FACS gates, one 16-nt barcode per
variant, Kd fitted per barcode from the mean bin (`compute_binding_Kd.Rmd` in each study's `extra/`).
The block's mean bin fit (`hill_gate_index`) is the matching method.

Their method, as the block reads it:
  - mean bin per barcode and concentration from cell-normalized counts, the sort fraction's role here.
    A gate with fewer reads than sorted cells keeps its raw counts; the block normalizes every gate;
  - Hill exponent fixed at 1, response in [2, 3] (2022) or [1.5, 3] (2020), baseline in [1, 1.5];
  - Kd bounded below at the lowest non-zero concentration / 100, above at the highest * 10 (2022) or
    * 100 (2020);
  - 2022 only: a mean bin is dropped when two non-neighbouring gates each hold > 40% of the cells
    (`calc.meanbin`), the rule the block applies as `constants.BIMODAL_SHARE`;
  - a point needs more than 2 cells; the block's floor is 3 reads;
  - a barcode needs enough cells at 80% (2022) or 60% (2020) of its concentrations, and 2022 drops
    fits with a normalized mean square residual far above the median. The block has neither filter.

The answer key is the authors' per-barcode `log10Ka` (= -log10 Kd), NA where their filters dropped it.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Library:
    study: str
    library: str
    answer_key: str
    # The answer key's library column value for this library.
    key_library: str
    # log10 Ka (= -log10 Kd) of the highest and the lowest non-zero concentration: the titrated range.
    neg_log10_kd_lower: float
    neg_log10_kd_upper: float


LIBRARIES = {
    "starr2020-lib1": Library("starr2020", "lib1", "answer-key/binding_Kds.csv", "lib1", 6.0, 13.0),
    "starr2022-pool1A": Library("starr2022", "pool1A", "answer-key/bc_binding.csv", "pool1A", 6.0, 13.0),
}

# The barcodes ending in these two bases, about one in eight, are fitted and compared; the rest are
# summed into one filler variant per sample, so every gate's read depth, and with it every mean bin,
# is the full library's.
SAMPLE_ENDINGS = ("AA", "CG")
FILLER = "OTHER_BARCODES"
