"""Phillips 2021 (CR6261) through the block, against the paper, pinned.

Four cases: two campaigns (H1: libraries a, b, c; H9: x, y, z), each with two settings, on the authors'
read counts with sort fractions and measured gate fluorescence bound.

  defaults  the settings a new block ships with (R² >= 0.8 and SE <= 1 gates on, plateau floor 0.667),
            on the metadata's concentration labels (10^-7 ... 10^-12 M)
  paper     the settings closest to Phillips' method: plateau floor 0.4, complete titrations only, both
            gates off, on the authors' fitting concentrations (labels shifted 0.05 decades lower)

Why the paper case turns the gates off: Phillips states R² >= 0.8 and SE <= 1 filters, but the answer key
(`20210323_<campaign>_avg_unf.csv`) is their unfiltered intermediate, and the gates move the block away from
it. Gates on vs off, otherwise the paper settings:
  H1  allClipped r 0.9919 vs 0.9941; Good n 1725 r 0.9911 vs n 1757 r 0.9924
  H9  allClipped r 0.9430 vs 0.9469; Good n 1270 r 0.9761 vs n 1290 r 0.9853
Median |d| differs by at most 0.0004 between the two in each pair. The concentration shift moves every
log10 K_D by exactly 0.05; under the default settings it brings H1's median |d| from 0.070 to 0.023.

Unlike Adams, the default R² >= 0.8 gate fails a minority of in-range replicate fits: 453 of 5,137 for
H1 and 808 of 3,954 for H9. The paper settings bring this down to 75 and 149.

Known gap, H9: 185 variants the authors fit at -log10 K_D 7.6 to 8.7 (replicate mean) are Above range
here. Their largest observed rise is 0.67 decades or less, below the block's 0.70-decade rise threshold
(LF_RISE_THRESHOLD_DECADES), so the block reports a bound where Phillips reports a weak, low-amplitude
binder. This keeps H9's allClipped within-0.5 share at 0.87.

Each case pins in pins/<case>.json:
  - the whole output, as a digest; a move fails with a diff of what changed
  - class counts, per variant and per replicate fit (the R² gate's toll shows in the latter)
  - agreement with the published replicate means: every variant (both sides clipped to the paper's
    7..12 range) and the Good variants whose published mean is inside it
  - the clones the authors also built and titrated one by one (isogenic K_D by flow cytometry)
"""

from __future__ import annotations

import dataclasses

import polars as pl
import pytest
import study_phillips2021 as study
from golden_support import UPDATE, agreement, check_digest, close, digest, snapshot, update_or_load
from paper_phillips2021 import GERMLINE, ISOGENIC, MATURE

pytestmark = pytest.mark.golden

PAPER = dataclasses.replace(
    study.PARAMS, plateau_floor=0.4, require_all_concentrations=True, r2_gate_min=None, kd_se_gate_max=None
)
# case -> (campaign, params, authors' shifted concentrations)
CASES = {
    "phillips2021-h1-defaults": ("h1", study.PARAMS, False),
    "phillips2021-h1-paper": ("h1", PAPER, True),
    "phillips2021-h9-defaults": ("h9", study.PARAMS, False),
    "phillips2021-h9-paper": ("h9", PAPER, True),
}
# How far agreement with the paper may fall before a test fails, beyond the pinned value.
R_DROP = 0.005
MEDIAN_RISE = 0.01
# The mature antibody against its isogenic K_D: the closest pin the paper offers to a true value.
MATURE_TOLERANCE = 0.3


def _current(case: str) -> tuple[dict, pl.DataFrame]:
    campaign, params, shifted = CASES[case]
    out = study.run_block(campaign, params, shifted)
    pc = out["per_clonotype"]
    ours = {k: (v, c) for k, v, c in pc.select("clonotypeKey", -pl.col("kd").log10(), "affinityClass").iter_rows()}
    pub = study.published(campaign)
    iso = study.isogenic(campaign)
    clones = {}
    for clone in ISOGENIC:
        if clone.name not in iso:
            continue
        v, cls = ours.get(clone.genotype, (None, "absent"))
        clones[clone.name] = {
            "genotype": clone.genotype,
            "log10Kd": None if v is None else round(-v, 6),
            "affinityClass": cls,
            "isogenic": round(iso[clone.name], 4),
            "paperTiteSeq": round(-pub[clone.genotype], 4),
        }
    germ_v, germ_cls = ours.get(GERMLINE, (None, "absent"))
    snap = snapshot(pc, "clonotypeKey")
    current = {
        "digest": digest(snap),
        "classes": {
            "variant": dict(sorted(pc["affinityClass"].value_counts().iter_rows())),
            "replicate": dict(sorted(out["per_replicate"]["affinityClass"].value_counts().iter_rows())),
            "replicateLowR2": out["per_replicate"].filter(pl.col("fitFailureReason") == "low_r2").height,
        },
        "agreement": {
            # The paper's convention: every variant with a value, both sides pinned to the 7..12 range.
            "allClipped": agreement(
                [(study.clip(v), study.clip(pub[k])) for k, (v, _) in ours.items() if v is not None and k in pub]
            ),
            # What a user filters on: Good values where the paper's mean is inside the range.
            "good": agreement(
                [(v, pub[k]) for k, (v, c) in ours.items() if c == "Good" and k in pub and study.in_range(pub[k])]
            ),
        },
        "clones": clones,
        "germline": {
            "negLog10Kd": None if germ_v is None else round(germ_v, 6),
            "affinityClass": germ_cls,
        },
    }
    return current, snap


@pytest.fixture(scope="module", params=list(CASES))
def case(request):
    current, snap = _current(request.param)
    return request.param, current, snap, update_or_load(request.param, current, snap)


def test_the_whole_output_matches_its_pin(case):
    # Any change to any variant's Kd,app, SE or class fails here, with a diff of what moved.
    name, current, snap, pins = case
    check_digest(name, snap, pins)


def test_class_counts_match_their_pin(case):
    name, current, _, pins = case
    assert current["classes"] == pins["classes"]


@pytest.mark.parametrize("subset", ["allClipped", "good"])
def test_agreement_with_the_paper_does_not_fall(case, subset):
    # The accuracy guard: a change may move numbers, but not further from the paper.
    name, current, _, pins = case
    now, pinned = current["agreement"][subset], pins["agreement"][subset]
    assert now["r"] >= pinned["r"] - R_DROP, f"{name} {subset}: r {pinned['r']} -> {now['r']}"
    assert now["medianAbsDiff"] <= pinned["medianAbsDiff"] + MEDIAN_RISE, (
        f"{name} {subset}: median |d| {pinned['medianAbsDiff']} -> {now['medianAbsDiff']}"
    )


def test_isogenic_clones_match_their_pin(case):
    # The variants the authors measured one by one; each moved clone is named with the paper's values.
    name, current, _, pins = case
    moved = [
        f"{clone} ({now['genotype']}): {pins['clones'][clone]['log10Kd']} ({pins['clones'][clone]['affinityClass']})"
        f" -> {now['log10Kd']} ({now['affinityClass']}); isogenic {now['isogenic']}, Tite-Seq {now['paperTiteSeq']}"
        for clone, now in current["clones"].items()
        if not close(now["log10Kd"], pins["clones"][clone]["log10Kd"])
        or now["affinityClass"] != pins["clones"][clone]["affinityClass"]
    ]
    assert set(current["clones"]) == set(pins["clones"])
    assert not moved, f"{name}:\n" + "\n".join(moved)


@pytest.mark.skipif(UPDATE, reason="pins are being rewritten")
def test_the_mature_antibody_lands_near_its_isogenic_kd(case):
    # 11111111111 carries all 11 library mutations: isogenic log10 K_D -9.50 (H1), -10.18 (H9).
    name, current, _, _ = case
    mature = current["clones"]["som"]
    assert mature["genotype"] == MATURE
    assert mature["affinityClass"] == "Good", f"{name}: the mature antibody rates {mature['affinityClass']}"
    assert abs(mature["log10Kd"] - mature["isogenic"]) < MATURE_TOLERANCE, (
        f"{name}: mature {mature['log10Kd']} vs isogenic {mature['isogenic']}"
    )


@pytest.mark.skipif(UPDATE, reason="pins are being rewritten")
def test_the_germline_does_not_bind(case):
    # The germline ancestor binds neither antigen detectably (the paper pins it at the 10^-7 M floor).
    name, current, _, _ = case
    assert current["germline"]["affinityClass"] == "Above range", f"{name}: germline {current['germline']}"
