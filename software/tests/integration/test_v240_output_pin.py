"""The v2.4.0 per-clonotype output, frozen.

"hill_gate_index" is what existing users run. This test pins its numbers on the committed corpus
so any change to the shared pipeline shows up here as a diff, not as a silent re-analysis.

Values were read off `main` @ ac1da53 (v2.4.0) and must not be regenerated to make a change pass.
Classes and reasons are exact; the three floats carry a 1e-6 relative tolerance for last-bit
platform noise in scipy's optimiser.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from constants import FitParams
from pipeline import run

CORPUS_DIR = Path(__file__).resolve().parent.parent / "data" / "corpus"
HILL = FitParams(equation="hill_gate_index")

# clonotypeKey -> (kd, hillCoefficient, r2, affinityClass, fitFailureReason, kdOutOfRange)
BIN_MODE = {
    "C0_ONLY": (None, None, None, "Failed", "insufficient_reads", None),
    "F_FLAT": (None, None, None, "Failed", "convergence_failure", None),
    "F_HOOK": (None, None, None, "Failed", "non_monotonic_signal", None),
    "F_INSUF_P": (None, None, None, "Failed", "insufficient_points", None),
    "F_INSUF_R": (None, None, None, "Failed", "insufficient_reads", None),
    "F_LOW_R2": (6.73364103816606e-14, 2.75075287505583, -2.3226727208225384e-09, "Failed", "low_r2", True),
    "G_HIGH": (9.735675760452655e-08, 1.293670843168412, 0.9700251556439752, "Good", None, False),
    "G_LOW": (1.078136104179333e-09, 0.9205933906222774, 0.9898017279662203, "Good", None, False),
    "G_MID": (1.2076379684272958e-08, 0.9965953054346497, 0.9905325322179203, "Good", None, False),
    "K_LOW": (7.10872899801838e-12, 0.5731588444591406, 0.96447738187534, "Good", None, True),
    "P_NOISY": (3.145850477127651e-09, 0.5860317483773869, 0.7420841201378108, "Partial", None, False),
    "P_N_HIGH": (1.0004407265071305e-08, 9.999999995033846, 0.9789435444936678, "Partial", None, False),
}

NO_BIN_MODE = {
    "NB_FILLER": (3.3691743385836444e-13, 5.0346830530391475, -4.409805853811122e-13, "Failed", "low_r2", True),
    "NB_GOOD": (3.002548924105987e-08, 0.8208794423987299, 0.9995920171324131, "Good", None, False),
    "NB_HOOK": (None, None, None, "Failed", "non_monotonic_signal", None),
    "NB_LOW_R2": (1.0916233883139839e-10, 9.999999999964755, 0.15685801409490896, "Failed", "low_r2", False),
}


@pytest.fixture(scope="module")
def bin_rows() -> dict[str, dict]:
    reads = pl.read_parquet(CORPUS_DIR / "reads_bin_mode.parquet")
    return {r["clonotypeKey"]: r for r in run(reads, params=HILL)["per_clonotype"].iter_rows(named=True)}


@pytest.fixture(scope="module")
def no_bin_rows() -> dict[str, dict]:
    reads = pl.read_parquet(CORPUS_DIR / "reads_no_bin_mode.parquet")
    return {r["clonotypeKey"]: r for r in run(reads, params=HILL)["per_clonotype"].iter_rows(named=True)}


def _assert_row(row: dict, expected: tuple) -> None:
    kd, hill, r2, affinity, reason, out_of_range = expected
    assert row["affinityClass"] == affinity
    assert row["fitFailureReason"] == reason
    assert row["kdOutOfRange"] is out_of_range
    for name, want in (("kd", kd), ("hillCoefficient", hill), ("r2", r2)):
        if want is None:
            assert row[name] is None, name
        else:
            assert row[name] == pytest.approx(want, rel=1e-6, abs=1e-12), name


@pytest.mark.parametrize("key,expected", sorted(BIN_MODE.items()))
def test_bin_mode_row_is_unchanged(bin_rows, key, expected):
    _assert_row(bin_rows[key], expected)


@pytest.mark.parametrize("key,expected", sorted(NO_BIN_MODE.items()))
def test_no_bin_mode_row_is_unchanged(no_bin_rows, key, expected):
    _assert_row(no_bin_rows[key], expected)


def test_the_corpus_carries_no_extra_clonotypes(bin_rows, no_bin_rows):
    assert set(bin_rows) == set(BIN_MODE)
    assert set(no_bin_rows) == set(NO_BIN_MODE)
