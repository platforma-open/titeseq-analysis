"""The v2.4.0 per-clonotype output, frozen.

"hill_gate_index" is what existing users run. This test pins its numbers on the committed corpus
so any change to the shared pipeline shows up here as a diff, not as a silent re-analysis.

Values were read off `main` @ ac1da53 (v2.4.0) and must not be regenerated to make a change pass.
Classes and reasons are exact. K_D and R-squared carry a 1e-6 relative tolerance for last-bit
platform noise in scipy's optimiser. The Hill exponent of a Failed row is not pinned at all, for
the reason given in `_assert_row`.

The other three frames (meanBin, fittedMeanBin, concentration value) are pinned whole against
frames v2.4.0 wrote on the same corpus, stored in tests/data/v240_reference/. Row order is part of
the pin, because the workflow deduplicates on the bytes.

One departure from v2.4.0 is intended and pinned on its own: a signal with no variance.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from constants import FitParams
from pipeline import run

CORPUS_DIR = Path(__file__).resolve().parent.parent / "data" / "corpus"
REFERENCE_DIR = Path(__file__).resolve().parent.parent / "data" / "v240_reference"
OTHER_FRAMES = ("mean_bin", "fitted_mean_bin", "concentration_value")
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
def bin_outputs() -> dict[str, pl.DataFrame]:
    return run(pl.read_parquet(CORPUS_DIR / "reads_bin_mode.parquet"), params=HILL)


@pytest.fixture(scope="module")
def no_bin_outputs() -> dict[str, pl.DataFrame]:
    return run(pl.read_parquet(CORPUS_DIR / "reads_no_bin_mode.parquet"), params=HILL)


@pytest.fixture(scope="module")
def bin_rows(bin_outputs) -> dict[str, dict]:
    return {r["clonotypeKey"]: r for r in bin_outputs["per_clonotype"].iter_rows(named=True)}


@pytest.fixture(scope="module")
def no_bin_rows(no_bin_outputs) -> dict[str, dict]:
    return {r["clonotypeKey"]: r for r in no_bin_outputs["per_clonotype"].iter_rows(named=True)}


def _assert_row(row: dict, expected: tuple) -> None:
    kd, hill, r2, affinity, reason, out_of_range = expected
    assert row["affinityClass"] == affinity
    assert row["fitFailureReason"] == reason
    assert row["kdOutOfRange"] is out_of_range
    for name, want in (("kd", kd), ("hillCoefficient", hill), ("r2", r2)):
        if want is None:
            assert row[name] is None, name
            continue
        # A rejected fit explains none of the signal, so its exponent is unconstrained and
        # scipy stops at a different value on each platform: NB_FILLER reads 5.03 on arm64
        # macOS and 6.19 on x86_64 linux. K_D and R-squared stay put on both. Pin the exponent
        # only where the block reports a usable fit.
        if name == "hillCoefficient" and affinity == "Failed":
            assert row[name] is not None, name
            continue
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


@pytest.mark.parametrize("frame", OTHER_FRAMES)
def test_bin_mode_frame_is_unchanged(bin_outputs, frame):
    expected = pl.read_parquet(REFERENCE_DIR / f"bin_mode_{frame}.parquet")
    assert_frame_equal(bin_outputs[frame], expected, check_exact=False, rel_tol=1e-6, abs_tol=1e-12)


@pytest.mark.parametrize("frame", OTHER_FRAMES)
def test_no_bin_mode_frame_is_unchanged(no_bin_outputs, frame):
    expected = pl.read_parquet(REFERENCE_DIR / f"no_bin_mode_{frame}.parquet")
    assert_frame_equal(no_bin_outputs[frame], expected, check_exact=False, rel_tol=1e-6, abs_tol=1e-12)


def test_a_flat_signal_is_failed_with_no_kd():
    """The intended departure from v2.4.0, which rated this Partial with a K_D of 28.39 M.

    Equal reads in every gate give a flat mean bin. The fit converges on the amplitude bound, R² is
    undefined (no variance), and the curve is an artefact: no K_D, exponent or R² is reported.
    """
    concs = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
    reads = pl.DataFrame(
        [
            {
                "clonotypeKey": "FLAT",
                "sampleId": f"s_c{i}_b{b}",
                "concentrationStr": str(c),
                "concentration": c,
                "bin": b,
                "reads": 25,
            }
            for i, c in enumerate(concs)
            for b in (1, 2, 3, 4)
        ]
    )
    out = run(reads, params=HILL)
    row = out["per_clonotype"].row(0, named=True)
    assert (row["affinityClass"], row["fitFailureReason"]) == ("Failed", "low_r2")
    assert row["kd"] is None
    assert row["hillCoefficient"] is None
    assert row["r2"] is None
    assert out["fitted_mean_bin"].height == 0
