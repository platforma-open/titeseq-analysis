"""Starr 2020 (lib1) and Starr 2022 (pool1A) through the block's mean bin fit, against the papers, pinned.

These are the mean bin fit's golden cases; Adams 2016 and Phillips 2021 cover the log fluorescence fit.
Each case runs the block's mean bin defaults with the sort fraction bound, on a fixed sample of the
library's barcodes (see `study_starr`). Starr 2022's analysis drops split mean bins as the block does;
Starr 2020's does not, so that case also shows what the rule costs on data it was not written for.

The block has neither of Starr's barcode filters (a minimum cell count at 80% of the concentrations, a
normalized-residual cut), so it reports values for barcodes the authors left NA. Agreement is measured
on the barcodes both report.

Each case pins in pins/<case>.json:
  - the whole output, as a digest; a move fails with a diff of what changed
  - class counts, and how many non-zero points have no mean bin: split ones, and ones with no read in
    a sorted gate
  - agreement with the authors' log10 Ka: every barcode both report (both sides clipped to the titrated
    range), and the Good barcodes whose published value is inside it
"""

from __future__ import annotations

import polars as pl
import pytest
import study_starr as study
from golden_support import UPDATE, agreement, check_digest, digest, snapshot, update_or_load
from paper_starr import FILLER, LIBRARIES

pytestmark = pytest.mark.golden

# How far agreement with the paper may fall before a test fails, beyond the pinned value.
R_DROP = 0.005
MEDIAN_RISE = 0.01


def _current(case: str) -> tuple[dict, pl.DataFrame]:
    out = study.run_block(case)
    pc = out["per_clonotype"].filter(pl.col("clonotypeKey").cast(pl.Utf8) != FILLER)
    ours = {
        k: (v, c)
        for k, v, c in pc.select(
            pl.col("clonotypeKey").cast(pl.Utf8), -pl.col("kd").log10(), "affinityClass"
        ).iter_rows()
    }
    pub = study.published(case)
    no_mean_bin = out["mean_bin"].filter(pl.col("meanBin").is_null() & (pl.col("concentrationStr") != "0.0")).height
    snap = snapshot(pc, "clonotypeKey")
    current = {
        "digest": digest(snap),
        "classes": {
            "variant": dict(sorted(pc["affinityClass"].value_counts().iter_rows())),
            "noMeanBinPoints": no_mean_bin,
        },
        "agreement": {
            "allClipped": agreement(
                [
                    (study.clip(case, v), study.clip(case, pub[k]))
                    for k, (v, _) in ours.items()
                    if v is not None and k in pub
                ]
            ),
            "good": agreement(
                [(v, pub[k]) for k, (v, c) in ours.items() if c == "Good" and k in pub and study.in_range(case, pub[k])]
            ),
        },
    }
    return current, snap


@pytest.fixture(scope="module", params=list(LIBRARIES))
def case(request):
    current, snap = _current(request.param)
    return request.param, current, snap, update_or_load(request.param, current, snap)


def test_the_whole_output_matches_its_pin(case):
    # Any change to any barcode's Kd,app, SE or class fails here, with a diff of what moved.
    name, current, snap, pins = case
    check_digest(name, snap, pins)


def test_class_counts_match_their_pin(case):
    _, current, _, pins = case
    assert current["classes"] == pins["classes"]


@pytest.mark.skipif(UPDATE, reason="pins are being rewritten")
@pytest.mark.parametrize("subset", ["allClipped", "good"])
def test_agreement_with_the_paper_does_not_fall(case, subset):
    # The accuracy guard: a change may move numbers, but not further from the paper.
    name, current, _, pins = case
    now, pinned = current["agreement"][subset], pins["agreement"][subset]
    assert now["r"] >= pinned["r"] - R_DROP, f"{name} {subset}: r {pinned['r']} -> {now['r']}"
    assert now["medianAbsDiff"] <= pinned["medianAbsDiff"] + MEDIAN_RISE, (
        f"{name} {subset}: median |d| {pinned['medianAbsDiff']} -> {now['medianAbsDiff']}"
    )
