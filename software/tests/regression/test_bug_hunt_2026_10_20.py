"""Reproductions from the 2026-10-20 bug hunt.

Each test asserts the correct behaviour and is marked strict xfail while the
defect is open. A test that starts passing flips to XPASS and fails the suite,
so the marker must be removed together with the fix.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import polars as pl
import pytest
from helpers import LINEAR_HILL

from classify import classify
from constants import FitParams
from hill_fit import fit_one_clonotype, weighted_r2
from io_layer import InputValidationError
from pipeline import run

CONCS_M = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]


def _bin_rows(
    concs: list[float],
    kds: list[float],
    rng: np.random.Generator,
    *,
    reads: int = 400,
    sigma: float = 0.4,
    prefix: str = "K",
    extra: Callable[[int, float, int], dict] | None = None,
    all_bins: bool = False,
) -> list[dict]:
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


def test_constant_signal_clonotype_is_never_good():
    # A clonotype whose observed signal is identical at every concentration (for
    # example all of its reads in one bin) has zero weighted variance, so
    # weighted_r2 returns NaN. A non-finite R² clears no threshold, so the fit is
    # Partial at best — never Good. The fitter refuses the flat signal before the
    # classifier ever sees it: a constant signal has no dynamic range, so
    # fit_one_clonotype returns converged=False and the pipeline reports Failed.
    y = np.full(6, 2.0)
    r2 = weighted_r2(y, y + 1e-3, np.full(6, 5.0))
    assert np.isnan(r2)
    assert classify(r2, 1.0, True, FitParams()).affinity_class == "Partial"

    flat = fit_one_clonotype(np.array(CONCS_M), y, np.full(6, 5.0), 2.0, bin_mode=True, max_bin_label=4, pin_n=True)
    assert not flat.converged
    assert flat.reason == "insufficient_dynamic_range"


def test_unbinned_zero_control_does_not_shift_kd():
    rng = np.random.default_rng(3)
    concs = [0.0, *CONCS_M]
    kds = list(10 ** rng.uniform(-9, -7, size=20))
    rows = _bin_rows(concs, kds, rng)
    ref = run(pl.DataFrame(rows), params=LINEAR_HILL)["per_clonotype"]
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
    out = run(unbinned, params=LINEAR_HILL)["per_clonotype"]
    joined = ref.select("clonotypeKey", pl.col("kd").alias("kd_ref")).join(
        out.select("clonotypeKey", "kd"), on="clonotypeKey"
    )
    ratio = joined.drop_nulls().select((pl.col("kd") / pl.col("kd_ref")).alias("r"))["r"]
    assert ratio.min() > 0.5  # ty: ignore[unsupported-operator]  polars types Series.min() as any scalar


def test_kd_not_capped_when_concentrations_are_nanomolar():
    rng = np.random.default_rng(5)
    concs_nm = [0.0, 0.3, 3.0, 30.0, 300.0, 3000.0, 10000.0]
    out = run(pl.DataFrame(_bin_rows(concs_nm, [5000.0] * 5, rng, reads=1000)), params=LINEAR_HILL)["per_clonotype"]
    kds = out["kd"].drop_nulls()
    assert kds.len() > 0
    assert (kds > 1500.0).all()


def test_sort_fraction_with_antigen_column_is_accepted():
    rng = np.random.default_rng(5)
    fr = {1: 0.4, 2: 0.3, 3: 0.2, 4: 0.1}
    rows = []
    for ag in ("AgA", "AgB"):
        part = _bin_rows(
            [0.0, *CONCS_M],  # the Adams model needs a 0 M control
            [1e-8] * 5,
            rng,
            all_bins=True,
            extra=lambda i, c, b, ag=ag: {"antigen": ag, "sort_fraction": fr[b]},
        )
        for r in part:
            r["sampleId"] = f"{ag}_{r['sampleId']}"
        rows += part
    run(
        pl.DataFrame(rows),
        params=FitParams(),
        target_antigen="AgA",
        antigen_column_ref="antigen",
        sort_fraction_col="sort_fraction",
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
        run(df, params=LINEAR_HILL, sort_fraction_col="sort_fraction")
    except InputValidationError as e:
        pytest.fail(f"rejected: {e}")


def test_too_few_points_for_parameter_count_is_reported_as_insufficient_points():
    rng = np.random.default_rng(5)
    concs = list(np.logspace(-10, -5, 4))  # 4 points, no 0 M control -> 4 free parameters
    out = run(
        pl.DataFrame(_bin_rows(concs, [1e-8] * 10, rng)),
        params=FitParams(fit_model="hill_free_n", min_concentration_points=3),
    )
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


def _sample_table(concs: list[float], bins: tuple[int, ...] = (1, 2, 3, 4)) -> pl.DataFrame:
    """The per-sample table the workflow passes as samples.tsv.

    Built from the sorting design, not from the reads, which is the whole point: the
    workflow reads it off the bin and concentration metadata columns, so a sample that
    collected no reads from any clonotype is still listed. `_bin_rows` names samples
    `s{conc index}_b{bin}`, so the design is the full (concentration x bin) product.
    """
    return pl.DataFrame(
        [
            {"sampleId": f"s{i}_b{b}", "bin": b, "concentration": c, "concentrationStr": str(c)}
            for i, c in enumerate(concs)
            for b in bins
        ]
    )


def test_grid_warning_silent_when_bins_simply_have_zero_reads(capsys):
    """bug-hunt M3: a clonotype with no reads in a bin is a true zero, not a gap.

    Every (bin, concentration) sample here was measured, so the check must stay quiet
    even though the reads table has no row for many (clonotype, bin) pairs — at 0 M the
    high bins are empty for the whole library.
    """
    rng = np.random.default_rng(21)
    concs = [0.0, *CONCS_M[:5]]
    rows = _bin_rows(concs, list(10 ** rng.uniform(-9, -7, 10)), rng, reads=300, sigma=0.35)
    run(pl.DataFrame(rows), params=LINEAR_HILL, samples=_sample_table(concs))
    assert "(R5)" not in capsys.readouterr().out


def test_grid_warning_fires_when_a_sample_is_missing_from_the_titration(capsys):
    """The R5 warning still reports a (bin, concentration) sample nobody measured."""
    rng = np.random.default_rng(21)
    concs = [0.0, *CONCS_M[:5]]
    rows = _bin_rows(concs, list(10 ** rng.uniform(-9, -7, 10)), rng, reads=300, sigma=0.35)
    samples = _sample_table(concs).filter(~((pl.col("bin") == 2) & (pl.col("concentration") == CONCS_M[3])))
    run(pl.DataFrame(rows), params=LINEAR_HILL, samples=samples)
    out = capsys.readouterr().out
    assert "(bin, concentration) sample(s) absent from the titration" in out
    assert "(R5)" in out
