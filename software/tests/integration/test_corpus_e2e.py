"""End-to-end tests against the committed deterministic corpus.

Load the reads parquets + manifest from tests/data/corpus and assert that the
full pipeline (io_layer → normalization → pre_fit → hill_fit → classify →
output_build) produces outputs matching each clonotype's declared expectation.

Corpus is regenerated manually via `tests/fixtures/generate_corpus.py`.
Run from blocks/titeseq-analysis/software/:

    uv sync
    uv run pytest tests/integration/test_corpus_e2e.py
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest

from pipeline import run

CORPUS_DIR = Path(__file__).resolve().parent.parent / "data" / "corpus"
MANIFEST = json.loads((CORPUS_DIR / "manifest.json").read_text())


def _entries(section: str) -> list[tuple[str, dict]]:
    """Flatten one corpus section's clonotype dict into (key, entry) tuples for parametrize."""
    return sorted(MANIFEST[section]["clonotypes"].items())


BIN_ENTRIES = _entries("bin_mode")
NO_BIN_ENTRIES = _entries("no_bin_mode")
ANTIGEN_ENTRIES = _entries("antigen")
ALL_ENTRIES = (
    [("bin_mode", k, e) for k, e in BIN_ENTRIES]
    + [("no_bin_mode", k, e) for k, e in NO_BIN_ENTRIES]
    + [("antigen", k, e) for k, e in ANTIGEN_ENTRIES]
)


@pytest.fixture(scope="module")
def bin_outputs() -> dict[str, pl.DataFrame]:
    reads = pl.read_parquet(CORPUS_DIR / "reads_bin_mode.parquet")
    return run(reads)


@pytest.fixture(scope="module")
def no_bin_outputs() -> dict[str, pl.DataFrame]:
    reads = pl.read_parquet(CORPUS_DIR / "reads_no_bin_mode.parquet")
    return run(reads)


@pytest.fixture(scope="module")
def antigen_outputs() -> dict[str, pl.DataFrame]:
    reads = pl.read_parquet(CORPUS_DIR / "reads_antigen.parquet")
    return run(
        reads,
        target_antigen=MANIFEST["antigen"]["target_antigen"],
        antigen_column_ref="antigen",
    )


@pytest.fixture(scope="module")
def outputs_by_section(
    bin_outputs, no_bin_outputs, antigen_outputs
) -> dict[str, dict[str, pl.DataFrame]]:
    return {
        "bin_mode": bin_outputs,
        "no_bin_mode": no_bin_outputs,
        "antigen": antigen_outputs,
    }


def _row_for(frame: pl.DataFrame, clonotype: str) -> dict:
    sub = frame.filter(pl.col("clonotypeKey") == clonotype)
    assert sub.height == 1, f"expected exactly one row for {clonotype}, found {sub.height}"
    return sub.row(0, named=True)


def _expected_classes(entry: dict) -> set[str]:
    if "expected_class" in entry:
        return {entry["expected_class"]}
    return set(entry["expected_class_in"])


def _expected_reasons(entry: dict) -> set:
    if "expected_reason_in" in entry:
        return set(entry["expected_reason_in"])
    return {entry["expected_reason"]}


# ---------------------------------------------------------------------------
# Parameterized behavioral checks — one pytest case per clonotype per section.
# When any of these fails, the pytest ID shows exactly which clonotype regressed,
# instead of a single giant failure listing all bad keys.
# ---------------------------------------------------------------------------


# Guards the R12 classification truth table: every declared clonotype lands in
# the expected class/reason branch regardless of corpus section.
@pytest.mark.parametrize(
    "section,key,entry",
    ALL_ENTRIES,
    ids=[f"{s}-{k}" for s, k, _ in ALL_ENTRIES],
)
def test_classification_matches_manifest(outputs_by_section, section, key, entry):
    row = _row_for(outputs_by_section[section]["per_clonotype"], key)
    assert row["affinityClass"] in _expected_classes(entry), (
        f"class={row['affinityClass']}, expected in {_expected_classes(entry)}"
    )
    assert row["fitFailureReason"] in _expected_reasons(entry), (
        f"reason={row['fitFailureReason']}, expected in {_expected_reasons(entry)}"
    )


# Tolerance bands guard against silent numeric drift in the Hill fit. Loose ranges
# let scipy version changes pass; a 10x kd shift would still fail.
@pytest.mark.parametrize(
    "section,key,entry",
    [(s, k, e) for s, k, e in ALL_ENTRIES if e.get("kd_range") is not None],
    ids=[
        f"{s}-{k}"
        for s, k, e in ALL_ENTRIES
        if e.get("kd_range") is not None
    ],
)
def test_kd_within_tolerance(outputs_by_section, section, key, entry):
    row = _row_for(outputs_by_section[section]["per_clonotype"], key)
    lo, hi = entry["kd_range"]
    assert row["kd"] is not None, "expected fitted kd, got null"
    assert lo <= row["kd"] <= hi, f"kd={row['kd']} not in [{lo}, {hi}]"


@pytest.mark.parametrize(
    "section,key,entry",
    [(s, k, e) for s, k, e in ALL_ENTRIES if e.get("hill_range") is not None],
    ids=[
        f"{s}-{k}"
        for s, k, e in ALL_ENTRIES
        if e.get("hill_range") is not None
    ],
)
def test_hill_within_tolerance(outputs_by_section, section, key, entry):
    row = _row_for(outputs_by_section[section]["per_clonotype"], key)
    lo, hi = entry["hill_range"]
    assert row["hillCoefficient"] is not None, "expected fitted n, got null"
    assert lo <= row["hillCoefficient"] <= hi, (
        f"n={row['hillCoefficient']} not in [{lo}, {hi}]"
    )


# R14b: kd outside [min non-zero conc, max conc] sets the flag. Skips entries
# that don't pin the flag to a specific value (kd_out_of_range == None in manifest).
@pytest.mark.parametrize(
    "section,key,entry",
    [(s, k, e) for s, k, e in ALL_ENTRIES if e.get("kd_out_of_range") is not None],
    ids=[
        f"{s}-{k}"
        for s, k, e in ALL_ENTRIES
        if e.get("kd_out_of_range") is not None
    ],
)
def test_kd_out_of_range_flag(outputs_by_section, section, key, entry):
    row = _row_for(outputs_by_section[section]["per_clonotype"], key)
    assert row["kdOutOfRange"] == entry["kd_out_of_range"]


# R17: null kd → kdPlotPosition = max_conc * 10; null n → hillPlotPosition = 1.0.
@pytest.mark.parametrize(
    "section,key,entry",
    [
        (s, k, e)
        for s, k, e in ALL_ENTRIES
        if e.get("hill_plot_position_is_sentinel")
    ],
    ids=[
        f"{s}-{k}"
        for s, k, e in ALL_ENTRIES
        if e.get("hill_plot_position_is_sentinel")
    ],
)
def test_r17_sentinel_values(outputs_by_section, section, key, entry):
    row = _row_for(outputs_by_section[section]["per_clonotype"], key)
    max_conc = MANIFEST["max_non_zero_concentration"]
    assert row["hillCoefficient"] is None
    assert row["kd"] is None
    assert row["hillPlotPosition"] == pytest.approx(1.0)
    assert row["kdPlotPosition"] == pytest.approx(max_conc * 10.0)


# ---------------------------------------------------------------------------
# Cross-section invariants — one-off assertions that aren't per-clonotype.
# ---------------------------------------------------------------------------


class TestFittedMeanBinMembership:
    """fitted_mean_bin must contain exactly the clonotypes with a successful fit."""

    @pytest.mark.parametrize(
        "section,fixture_name",
        [
            ("bin_mode", "bin_outputs"),
            ("no_bin_mode", "no_bin_outputs"),
        ],
    )
    def test_fitted_frame_membership(self, request, section, fixture_name):
        outputs = request.getfixturevalue(fixture_name)
        actual = set(outputs["fitted_mean_bin"]["clonotypeKey"].unique().to_list())
        expected = set(MANIFEST[section]["fitted_mean_bin_clonotypes"])
        assert actual == expected


class TestFrameInvariants:
    """Structural invariants on the output frames."""

    # fitted_mean_bin is the Hill-curve frame — defined only for c > 0.
    def test_fitted_mean_bin_excludes_c0(self, bin_outputs):
        fitted = bin_outputs["fitted_mean_bin"]
        assert fitted.filter(pl.col("concentration") == 0.0).height == 0

    def test_per_clonotype_row_count(self, bin_outputs):
        per = bin_outputs["per_clonotype"]
        expected_keys = set(MANIFEST["bin_mode"]["clonotypes"].keys())
        assert set(per["clonotypeKey"].to_list()) == expected_keys
        assert per.height == len(expected_keys)


class TestAntigenFilter:
    """R4: when target_antigen is specified, distractor clonotypes vanish from every frame."""

    # Same check across all three output frames — parametrize by frame name
    # so a leak into one frame doesn't mask leaks in the others.
    @pytest.mark.parametrize("frame_name", ["per_clonotype", "mean_bin", "fitted_mean_bin"])
    def test_distractors_filtered(self, antigen_outputs, frame_name):
        frame = antigen_outputs[frame_name]
        distractors = set(MANIFEST["antigen"]["distractor_clonotypes"])
        actual = set(frame["clonotypeKey"].to_list())
        assert actual.isdisjoint(distractors), (
            f"{frame_name}: distractors leaked: {actual & distractors}"
        )

    @pytest.mark.parametrize("frame_name", ["per_clonotype", "mean_bin", "fitted_mean_bin"])
    def test_targets_present(self, antigen_outputs, frame_name):
        frame = antigen_outputs[frame_name]
        targets = set(MANIFEST["antigen"]["clonotypes"].keys())
        actual = set(frame["clonotypeKey"].to_list())
        assert targets.issubset(actual), f"{frame_name}: missing targets {targets - actual}"


class TestDeterminism:
    """Re-running the pipeline must produce stable classifications + bounded
    numeric drift.

    Classification columns (affinityClass, fitFailureReason, kdOutOfRange) must
    match exactly. Fitted kd/n/r2 can drift slightly between runs because polars'
    multi-threaded group_by is not order-stable, which reshuffles the rows fed
    to scipy.optimize.curve_fit. For borderline fits (e.g. P_NOISY, which sits
    between Partial and Failed) the optimizer can follow a different path and
    land on a nearby local minimum, yielding a few percent drift in fitted
    values. Tolerance bands are set ~3x the empirically observed maxima so the
    test flags an order-of-magnitude regression in optimizer stability without
    flaking on routine FP-order noise.

    Empirical max drift across 15 run pairs: kd 0.01%, n 1.44%, r2 ~0.
    """

    def test_bin_mode_classification_repeatable(self):
        reads = pl.read_parquet(CORPUS_DIR / "reads_bin_mode.parquet")
        a = run(reads)["per_clonotype"].sort("clonotypeKey")
        b = run(reads)["per_clonotype"].sort("clonotypeKey")

        for col in ("clonotypeKey", "affinityClass", "fitFailureReason", "kdOutOfRange"):
            assert a[col].to_list() == b[col].to_list(), f"column {col} differs"

    def test_bin_mode_numeric_drift_bounded(self):
        reads = pl.read_parquet(CORPUS_DIR / "reads_bin_mode.parquet")
        a = run(reads)["per_clonotype"].sort("clonotypeKey")
        b = run(reads)["per_clonotype"].sort("clonotypeKey")
        keys = a["clonotypeKey"].to_list()

        # (column, rel, abs) — plot-position columns track their source values
        # when finite and exact sentinels when null, so the same bands apply.
        checks = [
            ("kd", 0.02, 1e-9),
            ("hillCoefficient", 0.05, 1e-9),
            ("r2", 0.0, 0.01),
            ("kdPlotPosition", 0.02, 1e-9),
            ("hillPlotPosition", 0.05, 1e-9),
        ]
        for col, rel, abs_tol in checks:
            for k, av, bv in zip(keys, a[col].to_list(), b[col].to_list()):
                if av is None and bv is None:
                    continue
                assert av is not None and bv is not None, (
                    f"{col} for {k}: one run null, other not ({av}, {bv})"
                )
                assert av == pytest.approx(bv, rel=rel, abs=abs_tol), (
                    f"{col} for {k}: {av} vs {bv} outside rel={rel}, abs={abs_tol}"
                )
