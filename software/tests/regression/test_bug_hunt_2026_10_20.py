"""Reproductions from the 2026-10-20 bug hunt.

Each test asserts the correct behaviour and is marked strict xfail while the
defect is open. A test that starts passing flips to XPASS and fails the suite,
so the marker must be removed together with the fix.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from classify import classify
from constants import FitParams
from hill_fit import weighted_r2
from io_layer import InputValidationError
from pipeline import run

CONCS_M = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]


def _bin_rows(concs, kds, rng, *, reads=400, sigma=0.4, prefix="K", extra=None, all_bins=False):
    rows = []
    for k, kd in enumerate(kds):
        for i, c in enumerate(concs):
            mb = 1.0 + 3.0 * (c / (kd + c) if c > 0 else 0.0)
            p = np.exp(-0.5 * ((np.arange(1, 5) - mb) / sigma) ** 2)
            p /= p.sum()
            cnt = rng.multinomial(reads, p)
            if all_bins:
                cnt = cnt + 1
            for b in range(1, 5):
                if cnt[b - 1] == 0:
                    continue
                row = {
                    "clonotypeKey": f"{prefix}{k:03d}",
                    "sampleId": f"s{i}_b{b}",
                    "concentration": float(c),
                    "bin": b,
                    "reads": int(cnt[b - 1]),
                }
                if extra:
                    row.update(extra(i, c, b))
                rows.append(row)
    return rows


@pytest.mark.xfail(strict=True, reason="bug-hunt H1: NaN R² falls through classify() to Partial")
def test_constant_signal_clonotype_is_not_partial():
    # A clonotype whose observed signal is identical at every concentration (for
    # example all of its reads in one bin) has zero weighted variance, so
    # weighted_r2 returns NaN. Whether the fit then converges is platform- and
    # float-noise-dependent; when it does, classify() must not pass it.
    y = np.full(6, 2.0)
    r2 = weighted_r2(y, y + 1e-3, np.full(6, 5.0))
    assert np.isnan(r2)
    assert classify(r2, 1.0, True, FitParams()).affinity_class == "Failed"


@pytest.mark.xfail(strict=True, reason="bug-hunt H2: unbinned 0 M control drags global baseline B to 0")
def test_unbinned_zero_control_does_not_shift_kd():
    rng = np.random.default_rng(3)
    concs = [0.0, *CONCS_M]
    kds = list(10 ** rng.uniform(-9, -7, size=20))
    rows = _bin_rows(concs, kds, rng)
    ref = run(pl.DataFrame(rows))["per_clonotype"]
    unbinned = (
        pl.DataFrame(rows)
        .with_columns(
            pl.when(pl.col("concentration") == 0)
            .then(pl.lit("s0_nobin"))
            .otherwise(pl.col("sampleId"))
            .alias("sampleId"),
            pl.when(pl.col("concentration") == 0).then(None).otherwise(pl.col("bin")).alias("bin"),
        )
        .group_by(["clonotypeKey", "sampleId", "concentration", "bin"], maintain_order=True)
        .agg(pl.col("reads").sum())
    )
    out = run(unbinned)["per_clonotype"]
    joined = ref.select("clonotypeKey", pl.col("kd").alias("kd_ref")).join(
        out.select("clonotypeKey", "kd"), on="clonotypeKey"
    )
    ratio = joined.drop_nulls().select((pl.col("kd") / pl.col("kd_ref")).alias("r"))["r"]
    assert ratio.min() > 0.5


@pytest.mark.xfail(strict=True, reason="bug-hunt H3: KD_HI = 1e3 caps Kd for nM concentrations")
def test_kd_not_capped_when_concentrations_are_nanomolar():
    rng = np.random.default_rng(5)
    concs_nm = [0.0, 0.3, 3.0, 30.0, 300.0, 3000.0, 10000.0]
    out = run(pl.DataFrame(_bin_rows(concs_nm, [5000.0] * 5, rng, reads=1000)))["per_clonotype"]
    kds = out["kd"].drop_nulls()
    assert kds.len() > 0
    assert (kds > 1500.0).all()


@pytest.mark.xfail(strict=True, reason="bug-hunt H4: sort-fraction sums checked before the antigen filter")
def test_sort_fraction_with_antigen_column_is_accepted():
    rng = np.random.default_rng(5)
    fr = {1: 0.4, 2: 0.3, 3: 0.2, 4: 0.1}
    rows = []
    for ag in ("AgA", "AgB"):
        part = _bin_rows(
            CONCS_M,
            [1e-8] * 5,
            rng,
            all_bins=True,
            extra=lambda i, c, b, ag=ag: {"antigen": ag, "sort_fraction": fr[b]},
        )
        for r in part:
            r["sampleId"] = f"{ag}_{r['sampleId']}"
        rows += part
    run(pl.DataFrame(rows), target_antigen="AgA", antigen_column_ref="antigen", sort_fraction_col="sort_fraction")


@pytest.mark.xfail(
    strict=True, reason="bug-hunt M1: a sample outside the titration with no sort fraction fails the run"
)
def test_sort_fraction_ignores_samples_without_concentration():
    rng = np.random.default_rng(5)
    rows = _bin_rows(CONCS_M, [1e-8] * 5, rng, all_bins=True, extra=lambda i, c, b: {"sort_fraction": 0.25})
    for k in range(5):
        rows.append(
            {
                "clonotypeKey": f"K{k:03d}",
                "sampleId": "library_presort",
                "concentration": None,
                "bin": None,
                "reads": 500,
                "sort_fraction": None,
            }
        )
    df = pl.DataFrame(
        rows, schema_overrides={"concentration": pl.Float64, "bin": pl.Int64, "sort_fraction": pl.Float64}
    )
    try:
        run(df, sort_fraction_col="sort_fraction")
    except InputValidationError as e:
        pytest.fail(f"rejected: {e}")


@pytest.mark.xfail(
    strict=True,
    reason="bug-hunt M2: 3-4 points pass the gate at minConcentrationPoints=3 and always fail as convergence_failure",
)
def test_too_few_points_for_parameter_count_is_reported_as_insufficient_points():
    rng = np.random.default_rng(5)
    concs = list(np.logspace(-10, -5, 4))  # 4 points, no 0 M control -> 4 free parameters
    out = run(pl.DataFrame(_bin_rows(concs, [1e-8] * 10, rng)), params=FitParams(min_concentration_points=3))
    reasons = set(out["per_clonotype"]["fitFailureReason"].to_list())
    assert "convergence_failure" not in reasons


@pytest.mark.xfail(
    strict=True, reason="bug-hunt H5: no-bin delta = 0.05 is unreachable for libraries above ~20 clonotypes"
)
def test_no_bin_mode_fits_binders_in_a_200_clonotype_library():
    rng = np.random.default_rng(9)
    n = 200
    abund = rng.lognormal(0, 0.5, n)
    abund /= abund.sum()
    kds = 10 ** rng.uniform(-9, -7, n)
    rows = []
    for i, c in enumerate([0.0, *CONCS_M]):
        p = abund * (c / (kds + c) + 0.02)
        p /= p.sum()
        cnt = rng.multinomial(200_000, p)
        for k in range(n):
            if cnt[k] > 0:
                rows.append(
                    {"clonotypeKey": f"K{k:03d}", "sampleId": f"s{i}", "concentration": c, "reads": int(cnt[k])}
                )
    pc = run(pl.DataFrame(rows))["per_clonotype"]
    assert (pc["affinityClass"] != "Failed").sum() > 0


@pytest.mark.xfail(strict=True, reason="bug-hunt M3: R5 grid warning fires on every sparse (MiXCR-shaped) dataset")
def test_grid_warning_silent_when_bins_simply_have_zero_reads(capsys):
    rng = np.random.default_rng(21)
    rows = _bin_rows([0.0, *CONCS_M[:5]], list(10 ** rng.uniform(-9, -7, 10)), rng, reads=300, sigma=0.35)
    run(pl.DataFrame(rows))
    assert "missing (bin, concentration) combinations" not in capsys.readouterr().out
