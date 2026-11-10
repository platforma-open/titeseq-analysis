"""The "log_fluorescence" per-clonotype output, frozen.

"log_fluorescence" is the equation a new project runs. This test pins its numbers on the committed
corpus, so a change to the fit, the gates or the classes shows up here as a diff, not as a silent
re-analysis. It needs no study data and runs in every suite, CI included.

Two settings are pinned: the block's defaults for a new project (model/src/fitSettings.ts,
kind/src/index.ts), and the same with both fit-quality gates off, which keeps the classes the gates
override in view. The equation needs gates, so only the bin-mode corpus is run. The corpus has no
synonymous variants, so there is no synonymous output to pin.

Values are the block's own output on the committed corpus, recorded when this pin was added, and must
not be regenerated to make a change pass. Classes and reasons are exact. K_D, its SE and R² carry a
1e-6 relative tolerance for last-bit platform noise in scipy's optimiser.

The other three frames (meanBin, fittedMeanBin, concentration value) are pinned whole against
frames written by the same code, stored in tests/data/log_fluorescence_reference/. Row order is
part of the pin, because the workflow deduplicates on the bytes.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from constants import FitParams
from pipeline import run

CORPUS_DIR = Path(__file__).resolve().parent.parent / "data" / "corpus"
REFERENCE_DIR = Path(__file__).resolve().parent.parent / "data" / "log_fluorescence_reference"
OTHER_FRAMES = ("mean_bin", "fitted_mean_bin", "concentration_value")
PARAMS = {
    "defaults": FitParams(equation="log_fluorescence", plateau_floor=0.667, r2_gate_min=0.8, kd_se_gate_max=1.0),
    "gates_off": FitParams(equation="log_fluorescence", plateau_floor=0.667),
}

# clonotypeKey -> (kd, kdLog10Se, r2, affinityClass, fitFailureReason, kdOutOfRange)
DEFAULTS = {
    "C0_ONLY": (None, None, None, "Failed", "insufficient_reads", None),
    "F_FLAT": (1e-06, None, None, "Above range", None, True),
    "F_HOOK": (1.4216854477623837e-08, 0.6026460191898025, 0.6542655556747916, "Failed", "low_r2", False),
    "F_INSUF_P": (None, None, None, "Failed", "insufficient_points", None),
    "F_INSUF_R": (None, None, None, "Failed", "insufficient_reads", None),
    "F_LOW_R2": (1e-10, None, None, "Below range", None, True),
    "G_HIGH": (7.217910926273202e-07, 0.5103237629535269, 0.9042278488613067, "Partial", None, False),
    "G_LOW": (5.798103553042284e-09, 0.08313580129537443, 0.9841121556494522, "Good", None, False),
    "G_MID": (7.495039599187662e-08, 0.1542689340801187, 0.9750505637670569, "Good", None, False),
    "K_LOW": (1e-10, None, None, "Below range", None, True),
    "P_NOISY": (1.245993582351643e-08, 0.3408202225376994, 0.7701473890417541, "Failed", "low_r2", False),
    "P_N_HIGH": (9.870984131013337e-08, 0.46535617831420345, 0.8613249743873173, "Good", None, False),
}
# The gates only rate a fit; with them off, the two rows they failed keep the fit's own class.
GATES_OFF = {
    **DEFAULTS,
    "F_HOOK": (1.4216854477623837e-08, 0.6026460191898025, 0.6542655556747916, "Partial", None, False),
    "P_NOISY": (1.245993582351643e-08, 0.3408202225376994, 0.7701473890417541, "Good", None, False),
}
EXPECTED = {"defaults": DEFAULTS, "gates_off": GATES_OFF}


@pytest.fixture(scope="module")
def reads() -> pl.DataFrame:
    return pl.read_parquet(CORPUS_DIR / "reads_bin_mode.parquet")


@pytest.fixture(scope="module", params=list(PARAMS))
def case(request, reads) -> tuple[str, dict[str, pl.DataFrame]]:
    return request.param, run(reads, params=PARAMS[request.param])


def test_each_row_is_unchanged(case):
    name, outputs = case
    rows = {r["clonotypeKey"]: r for r in outputs["per_clonotype"].iter_rows(named=True)}
    assert set(rows) == set(EXPECTED[name])
    moved = []
    for key, (kd, se, r2, affinity, reason, out_of_range) in sorted(EXPECTED[name].items()):
        row = rows[key]
        if (row["affinityClass"], row["fitFailureReason"], row["kdOutOfRange"]) != (affinity, reason, out_of_range):
            moved.append(f"{key}: {row['affinityClass']}/{row['fitFailureReason']}/{row['kdOutOfRange']}")
        for col, want in (("kd", kd), ("kdLog10Se", se), ("r2", r2)):
            got = row[col]
            if want is None or got is None:
                if got is not want:
                    moved.append(f"{key} {col}: {want} -> {got}")
            elif got != pytest.approx(want, rel=1e-6, abs=1e-12):
                moved.append(f"{key} {col}: {want} -> {got}")
    assert not moved, f"{name}:\n" + "\n".join(moved)


def test_the_equation_fixes_the_exponent_at_one(case):
    _, outputs = case
    fitted = outputs["per_clonotype"].filter(pl.col("kd").is_not_null())
    assert fitted["hillCoefficient"].to_list() == [1.0] * fitted.height


@pytest.mark.parametrize("frame", OTHER_FRAMES)
def test_frame_is_unchanged(case, frame):
    name, outputs = case
    expected = pl.read_parquet(REFERENCE_DIR / f"{name}_{frame}.parquet")
    assert_frame_equal(outputs[frame], expected, check_exact=False, rel_tol=1e-6, abs_tol=1e-12)
