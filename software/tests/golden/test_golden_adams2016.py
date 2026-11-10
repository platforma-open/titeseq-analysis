"""Adams 2016 through the block, against the paper, pinned.

The authors' read counts, 3 replicates, sort fractions, synonymous mode, in two cases:
  paper     the settings closest to the paper: the block's defaults with both fit-quality gates off,
            as the authors filter no fit
  defaults  the settings a new block ships with (R² >= 0.8 and SE <= 1 gates on, plateau floor 0.667),
            what a user gets without changing a setting

Each case pins three things in pins/<case>.json:
  - the whole output, as a digest; a move fails with a diff of what changed
  - agreement with the authors' values (their per-replicate fits, clipped and averaged as their code does)
  - the Table I clones, the variants the paper also measured by flow cytometry
"""

from __future__ import annotations

import dataclasses

import polars as pl
import pytest
import study_adams2016 as study
from golden_support import UPDATE, agreement, check_digest, close, digest, snapshot, update_or_load
from paper_adams2016 import LOG10_KD_LOWER, LOG10_KD_UPPER, TABLE_I

pytestmark = pytest.mark.golden

CASES = {
    "adams2016-paper": dataclasses.replace(study.PARAMS, r2_gate_min=None, kd_se_gate_max=None),
    "adams2016-defaults": study.PARAMS,
}
# How far agreement with the paper may fall before a test fails, beyond the pinned value.
R_DROP = 0.005
MEDIAN_RISE = 0.01


def _log10(frame: pl.DataFrame, key: str) -> dict[str, tuple[float | None, str]]:
    return {
        k: (v, c)
        for k, v, c in frame.select(
            pl.col(key).cast(pl.Utf8), pl.col("kd").log10(), pl.col("affinityClass")
        ).iter_rows()
    }


def _clip(x: float) -> float:
    return min(max(x, LOG10_KD_LOWER), LOG10_KD_UPPER)


def _current(case: str) -> tuple[dict, pl.DataFrame]:
    out = study.run_block(CASES[case])
    nt = _log10(out["per_clonotype"], "clonotypeKey")
    aa = _log10(out["per_variant"], "aaKey")
    pub_nt, pub_aa = study.published()

    def metrics(ours: dict, pub: dict) -> dict:
        return {
            # The paper's convention: every variant with a value, both sides clipped to the limits.
            "allClipped": agreement([(_clip(v), pub[k]) for k, (v, _) in ours.items() if v is not None and k in pub]),
            # What a user filters on: Good values where the paper measured inside the range.
            "good": agreement(
                [(v, pub[k]) for k, (v, c) in ours.items() if c == "Good" and k in pub and study.in_range(pub[k])]
            ),
        }

    keys = study.clone_keys()
    clones = {}
    for c in TABLE_I:
        v, cls = nt.get(keys.get(c.name, ""), (None, "absent"))
        clones[c.name] = {
            "log10Kd": None if v is None else round(v, 6),
            "affinityClass": cls,
            "paperTiteSeq": c.tite_seq_log10_kd,
            "paperFlow": c.flow_log10_kd,
        }
    snap = pl.concat(
        [
            snapshot(out["per_clonotype"], "clonotypeKey").with_columns(("nt:" + pl.col("key")).alias("key")),
            snapshot(out["per_variant"], "aaKey").with_columns(("aa:" + pl.col("key")).alias("key")),
        ]
    ).sort("key")
    current = {
        "digest": digest(snap),
        "classes": {
            "nt": dict(sorted(out["per_clonotype"]["affinityClass"].value_counts().iter_rows())),
            "aa": dict(sorted(out["per_variant"]["affinityClass"].value_counts().iter_rows())),
        },
        "agreement": {"nt": metrics(nt, pub_nt), "aa": metrics(aa, pub_aa)},
        "clones": clones,
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


@pytest.mark.parametrize("grain", ["nt", "aa"])
@pytest.mark.parametrize("subset", ["allClipped", "good"])
def test_agreement_with_the_paper_does_not_fall(case, grain, subset):
    # The accuracy guard: a change may move numbers, but not further from the paper.
    name, current, _, pins = case
    now, pinned = current["agreement"][grain][subset], pins["agreement"][grain][subset]
    assert now["r"] >= pinned["r"] - R_DROP, f"{name} {grain} {subset}: r {pinned['r']} -> {now['r']}"
    assert now["medianAbsDiff"] <= pinned["medianAbsDiff"] + MEDIAN_RISE, (
        f"{name} {grain} {subset}: median |d| {pinned['medianAbsDiff']} -> {now['medianAbsDiff']}"
    )


def test_table_i_clones_match_their_pin(case):
    # The variants the paper measured one by one; each moved clone is named with the paper's values.
    name, current, _, pins = case
    moved = [
        f"{clone}: {pins['clones'][clone]['log10Kd']} ({pins['clones'][clone]['affinityClass']}) -> "
        f"{now['log10Kd']} ({now['affinityClass']}); paper Tite-Seq {now['paperTiteSeq']}, flow {now['paperFlow']}"
        for clone, now in current["clones"].items()
        if not close(now["log10Kd"], pins["clones"][clone]["log10Kd"])
        or now["affinityClass"] != pins["clones"][clone]["affinityClass"]
    ]
    assert not moved, f"{name}:\n" + "\n".join(moved)


@pytest.mark.skipif(UPDATE, reason="pins are being rewritten")
def test_the_wild_type_lands_near_the_paper(case, request):
    # pJK36 is the wild type: Adams reports -8.92 by Tite-Seq, -8.61 by flow cytometry.
    name, current, _, _ = case
    if name == "adams2016-defaults":
        request.applymarker(
            pytest.mark.xfail(strict=True, reason="the shipped R² >= 0.8 gate fails the wild type on Adams' counts")
        )
    wt = current["clones"]["pJK36"]
    assert wt["log10Kd"] is not None, f"{name}: the wild type has no Kd,app"
    assert abs(wt["log10Kd"] - wt["paperTiteSeq"]) < 0.5
