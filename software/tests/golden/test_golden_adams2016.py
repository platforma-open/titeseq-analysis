"""Adams 2016 through the block, against the paper, pinned.

One case: the authors' read counts, 3 replicates, sort fractions, synonymous mode, with the settings
closest to the paper: the block's defaults with both fit-quality gates off, as the authors filter no fit.
(The shipped defaults are pinned on Phillips 2021. On Adams' counts their R² >= 0.8 gate fails most
in-range fits, the wild type's included.)

It pins four things in pins/<case>.json:
  - the whole output, as a digest; a move fails with a diff of what changed. It covers the
    per-nucleotide-variant table (the outlier test's columns) and each amino-acid variant's outlier note
  - agreement with the authors' values (their per-replicate fits, clipped and averaged as their code does)
  - the Table I clones, the variants the paper also measured by flow cytometry
  - the synonymous outlier test: the dataset sigma, the counts, the flagged variants and the spread of z

A nucleotide variant key is its sequence, so the pins name a flagged variant by a hash of its key.
"""

from __future__ import annotations

import dataclasses
import hashlib
import statistics

import polars as pl
import pytest
import study_adams2016 as study
from golden_support import UPDATE, agreement, check_digest, close, digest, snapshot, update_or_load
from paper_adams2016 import LOG10_KD_LOWER, LOG10_KD_UPPER, TABLE_I

from constants import SYNONYM_OUTLIER_Z

pytestmark = pytest.mark.golden

CASES = {
    "adams2016-paper": dataclasses.replace(study.PARAMS, r2_gate_min=None, kd_se_gate_max=None),
}
# How far agreement with the paper may fall before a test fails, beyond the pinned value.
R_DROP = 0.005
MEDIAN_RISE = 0.01

# The SD of a standard normal kept inside |z| < SYNONYM_OUTLIER_Z: what the z of the in-range,
# non-flagged variants should spread by when sigma and the SEs are right.
_c = SYNONYM_OUTLIER_Z
_n = statistics.NormalDist()
Z_SD_NOT_FLAGGED = (1 - 2 * _c * _n.pdf(_c) / (2 * _n.cdf(_c) - 1)) ** 0.5
# 3 times the SD's sampling error on Adams, by a bootstrap over amino-acid variants (0.027).
Z_SD_TOLERANCE = 0.08


def _log10(frame: pl.DataFrame, key: str) -> dict[str, tuple[float | None, str]]:
    return {
        k: (v, c)
        for k, v, c in frame.select(
            pl.col(key).cast(pl.Utf8), pl.col("kd").log10(), pl.col("affinityClass")
        ).iter_rows()
    }


def _clip(x: float) -> float:
    return min(max(x, LOG10_KD_LOWER), LOG10_KD_UPPER)


def _key_hash(key: str) -> str:
    return hashlib.sha256(key.encode()).hexdigest()[:12]


def _snapshot(out: dict) -> pl.DataFrame:
    """Both grains' values, the per-nucleotide-variant table and the amino-acid outlier note."""
    per_nt = out["per_nt_variant"].with_columns(pl.col("ntKey").cast(pl.Utf8).alias("key")).drop("ntKey")
    per_nt = per_nt.with_columns(pl.col(c).round(6) for c, t in per_nt.schema.items() if t == pl.Float64)
    nt = snapshot(out["per_clonotype"], "clonotypeKey").join(per_nt, on="key", how="left")
    aa = snapshot(out["per_variant"], "aaKey").join(
        out["per_variant"].select(pl.col("aaKey").cast(pl.Utf8).alias("key"), "synonymOutliers"), on="key", how="left"
    )
    return pl.concat(
        [nt.with_columns(("nt:" + pl.col("key")).alias("key")), aa.with_columns(("aa:" + pl.col("key")).alias("key"))],
        how="diagonal",
    ).sort("key")


def _outliers(out: dict) -> dict:
    """The outlier test's summary. A variant at a wall is tested as a bound and kept out of the z spread."""
    summary = out["synonymous_summary"]["synonymOutlier"]
    rows = out["per_nt_variant"].join(
        out["per_clonotype"].select(pl.col("clonotypeKey").cast(pl.Utf8).alias("ntKey"), "kdOutOfRange"), on="ntKey"
    )
    tested = rows.filter(pl.col("kdLog10SynonymDelta").is_not_null())
    in_range = tested.filter(~pl.col("kdOutOfRange"))
    z = in_range["kdLog10SynonymZ"].to_list()
    z_not_flagged = [v for v in z if abs(v) < SYNONYM_OUTLIER_Z]
    flagged = tested.filter(pl.col("synonymOutlier"))
    return {
        "datasetSigma": round(summary["datasetSigma"], 6),
        "tested": {
            "nt": tested.height,
            "inRange": in_range.height,
            "atWall": tested.height - in_range.height,
            "aa": summary["testedAaVariants"],
        },
        "flagged": {
            "nt": flagged.height,
            "inRange": flagged.filter(~pl.col("kdOutOfRange")).height,
            "atWall": flagged.filter(pl.col("kdOutOfRange")).height,
            "aa": summary["flaggedAaVariants"],
        },
        "zSd": {
            "inRange": round(statistics.stdev(z), 6),
            "inRangeNotFlagged": round(statistics.stdev(z_not_flagged), 6),
        },
        "flaggedVariants": {
            _key_hash(key): {"delta": round(delta, 6), "z": round(zv, 6), "atWall": wall}
            for key, delta, zv, wall in flagged.select(
                "ntKey", "kdLog10SynonymDelta", "kdLog10SynonymZ", "kdOutOfRange"
            ).iter_rows()
        },
    }


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
        "synonymOutlier": _outliers(out),
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
def test_the_wild_type_lands_near_the_paper(case):
    # pJK36 is the wild type: Adams reports -8.92 by Tite-Seq, -8.61 by flow cytometry.
    name, current, _, _ = case
    wt = current["clones"]["pJK36"]
    assert wt["log10Kd"] is not None, f"{name}: the wild type has no Kd,app"
    assert abs(wt["log10Kd"] - wt["paperTiteSeq"]) < 0.5


def test_the_outlier_test_summary_matches_its_pin(case):
    name, current, _, pins = case
    now, pinned = current["synonymOutlier"], pins["synonymOutlier"]
    assert close(now["datasetSigma"], pinned["datasetSigma"]), (
        f"{name}: dataset sigma {pinned['datasetSigma']} -> {now['datasetSigma']}"
    )
    assert (now["tested"], now["flagged"]) == (pinned["tested"], pinned["flagged"])
    for subset, value in now["zSd"].items():
        assert close(value, pinned["zSd"][subset]), f"{name}: SD of z ({subset}) {pinned['zSd'][subset]} -> {value}"


def test_the_flagged_variants_match_their_pin(case):
    # Each moved, new or dropped flag is named by its key's hash.
    name, current, _, pins = case
    now, pinned = current["synonymOutlier"]["flaggedVariants"], pins["synonymOutlier"]["flaggedVariants"]
    moved = [f"{k}: newly flagged {now[k]}" for k in sorted(set(now) - set(pinned))]
    moved += [f"{k}: no longer flagged, was {pinned[k]}" for k in sorted(set(pinned) - set(now))]
    moved += [
        f"{k}: {pinned[k]} -> {now[k]}"
        for k in sorted(set(now) & set(pinned))
        if not close(now[k]["delta"], pinned[k]["delta"])
        or not close(now[k]["z"], pinned[k]["z"])
        or now[k]["atWall"] != pinned[k]["atWall"]
    ]
    assert not moved, f"{name}:\n" + "\n".join(moved)


def test_z_of_the_in_range_variants_spreads_as_a_standard_normal(case):
    # With sigma and the SEs right, z over the in-range variants not flagged is a standard normal
    # cut at the threshold. Too narrow: the errors are overstated, and real outliers go unflagged.
    # Too wide: they are understated, and agreeing variants get flagged.
    name, current, _, _ = case
    sd = current["synonymOutlier"]["zSd"]["inRangeNotFlagged"]
    assert abs(sd - Z_SD_NOT_FLAGGED) <= Z_SD_TOLERANCE, (
        f"{name}: SD of z over in-range, non-flagged variants {sd:.3f}, expected {Z_SD_NOT_FLAGGED:.3f}"
    )
