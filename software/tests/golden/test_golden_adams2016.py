"""Adams 2016 through the block, against the paper, pinned.

The authors' read counts, 3 replicates, sort fractions, each replicate's measured gate fluorescence,
synonymous mode, in two cases:
  paper     the settings closest to the paper: the block's defaults with both fit-quality gates off,
            as the authors filter no fit
  defaults  the settings a new block ships with (R² >= 0.8 and SE <= 1 gates on, plateau floor 0.667),
            what a user gets without changing a setting

A third run, `adams2016-parent-delta`, is the defaults with every amino-acid variant linked to the
wild type. It pins the parent delta of each variant and how many deltas were taken inside the
replicates.

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


def _snapshot(out: dict) -> pl.DataFrame:
    """Both grains' values and the per-nucleotide-variant table."""
    per_nt = out["per_nt_variant"].with_columns(pl.col("ntKey").cast(pl.Utf8).alias("key")).drop("ntKey")
    per_nt = per_nt.with_columns(pl.col(c).round(6) for c, t in per_nt.schema.items() if t == pl.Float64)
    nt = snapshot(out["per_clonotype"], "clonotypeKey").join(per_nt, on="key", how="left")
    aa = snapshot(out["per_variant"], "aaKey")
    return pl.concat(
        [nt.with_columns(("nt:" + pl.col("key")).alias("key")), aa.with_columns(("aa:" + pl.col("key")).alias("key"))],
        how="diagonal",
    ).sort("key")


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
    snap = _snapshot(out)
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
def test_agreement_with_the_authors_fits_does_not_fall(case, grain, subset):
    # The accuracy guard: a change may move numbers, but not further from the authors' fits. The
    # nucleotide key is their published average; the amino-acid key is built from their fits here.
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
def test_the_wild_type_lands_near_the_paper(case):
    # pJK36 is the wild type: Adams reports -8.92 by Tite-Seq, -8.61 by flow cytometry. Under the shipped
    # R² gate it passes because R² counts its 0 M point (`log_fluorescence.r2_with_zero_point`).
    name, current, _, _ = case
    wt = current["clones"]["pJK36"]
    assert wt["log10Kd"] is not None, f"{name}: the wild type has no Kd,app"
    assert abs(wt["log10Kd"] - wt["paperTiteSeq"]) < 0.5


# ---- The parent delta, against the wild type.

PARENT_CASE = "adams2016-parent-delta"


@pytest.fixture(scope="module")
def parent_delta():
    out = study.run_block(study.PARAMS, parent=True)
    # The snapshot's columns hold the delta and its SE, so a move is reported as the other pins are.
    snap = (
        out["per_variant"]
        .select(
            pl.col("aaKey").cast(pl.Utf8).alias("key"),
            pl.col("kdLog10ParentDelta").round(6).alias("log10Kd"),
            pl.col("kdLog10ParentDeltaSe").round(6).alias("log10KdSe"),
            pl.col("affinityClass"),
        )
        .sort("key")
    )
    summary = out["synonymous_summary"]
    current = {
        "digest": digest(snap),
        "replicateDeltas": summary["replicateDeltas"],
        "combinedValueDeltas": summary["combinedValueDeltas"],
        "combinedValueDeltasInRange": summary["combinedValueDeltasInRange"],
        "deltas": snap["log10Kd"].count(),
        "deltasWithSe": snap["log10KdSe"].count(),
    }
    return current, snap, update_or_load(PARENT_CASE, current, snap)


def test_every_parent_delta_matches_its_pin(parent_delta):
    _, snap, pins = parent_delta
    check_digest(PARENT_CASE, snap, pins)


def test_the_count_of_deltas_taken_inside_the_replicates_matches_its_pin(parent_delta):
    current, _, pins = parent_delta
    counted = {k: v for k, v in current.items() if k != "digest"}
    assert counted == {k: v for k, v in pins.items() if k != "digest"}


def test_the_parent_links_leave_every_kd_as_it_was(parent_delta):
    # The parent delta is read off the fits. It changes none of them.
    _, snap, _ = parent_delta
    linked = study.run_block(study.PARAMS, parent=True)["per_variant"]
    plain = study.run_block(study.PARAMS)["per_variant"]
    assert snapshot(linked, "aaKey").equals(snapshot(plain, "aaKey"))
