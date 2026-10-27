"""Bin-mode pipeline for the poisson_log model (Adams, Mora, Walczak, Kinney 2016, Appendix E).

Turns the long reads frame into count tensors, fits every clonotype with `poisson_fit`, and returns
the per-clonotype and fitted-mean-bin frames in the same schemas as the hill_free_n path.
"""

from __future__ import annotations

import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
import polars as pl

from constants import COL_BIN, COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL, COL_READS, COL_SAMPLE, FitParams
from cpus import available_cpus
from hill_fit import weighted_r2
from io_layer import InputValidationError
from log import log
from normalization import CLONOTYPE_READS_AT_CONC, SIGNAL
from output_build import FITTED_MEAN_BIN_SCHEMA, PER_CLONOTYPE_SCHEMA
from poisson_fit import (
    BASELINE_CONTAMINATION_RATIO,
    BIN_FLUORESCENCE_MAX_TOP_SPREAD,
    GridFit,
    baseline,
    exposure,
    fit_titration,
    predicted_mean_bin,
    predicted_mean_bin_measured,
)

_PARALLEL_MIN_CLONOTYPES = 400
_CHUNK = 200


@dataclass(frozen=True)
class CountTensors:
    clonotypes: list[str]
    conc_strs: list[str]
    concentrations: np.ndarray  # (C,) ascending, may include 0
    bin_labels: np.ndarray  # (B,) ascending
    reads: np.ndarray  # (V, C, B) raw reads
    sample_cells: np.ndarray | None  # (C, B) cells sorted, or None
    # (C, B) measured mean log10 fluorescence per sample, or None. NaN where a sample had no value.
    bin_log_fluorescence: np.ndarray | None


def build_tensors(
    reads: pl.DataFrame, cell_count_col: str | None, bin_fluorescence_col: str | None = None
) -> CountTensors:
    clonotypes = sorted(reads[COL_CLONOTYPE].unique().to_list())
    # Rows without a bin (an unsorted 0 M control) carry no bin information; the model counts sorted bins only.
    reads = reads.filter(pl.col(COL_BIN).is_not_null())
    concs = reads.select(COL_CONC_STR, COL_CONC_VAL).unique().sort(COL_CONC_VAL)
    bins = np.array(sorted(reads[COL_BIN].unique().to_list()))
    ci = {c: i for i, c in enumerate(concs[COL_CONC_STR].to_list())}
    bi = {int(b): i for i, b in enumerate(bins)}
    vi = {c: i for i, c in enumerate(clonotypes)}

    agg = reads.group_by([COL_CLONOTYPE, COL_CONC_STR, COL_BIN]).agg(pl.col(COL_READS).sum())
    R = np.zeros((len(clonotypes), concs.height, bins.size))
    np.add.at(
        R,
        (
            np.array([vi[x] for x in agg[COL_CLONOTYPE]]),
            np.array([ci[x] for x in agg[COL_CONC_STR]]),
            np.array([bi[int(x)] for x in agg[COL_BIN]]),
        ),
        agg[COL_READS].cast(pl.Float64).to_numpy(),
    )
    cells = None
    if cell_count_col is not None:
        per_sample = reads.select(COL_SAMPLE, COL_CONC_STR, COL_BIN, cell_count_col).unique()
        cells = np.zeros((concs.height, bins.size))
        for row in per_sample.iter_rows(named=True):
            cells[ci[row[COL_CONC_STR]], bi[int(row[COL_BIN])]] += float(row[cell_count_col])
    blf = None
    if bin_fluorescence_col is not None:
        per_sample = reads.select(COL_SAMPLE, COL_CONC_STR, COL_BIN, bin_fluorescence_col).unique()
        blf = np.full((concs.height, bins.size), np.nan)
        for row in per_sample.iter_rows(named=True):
            v = row[bin_fluorescence_col]
            if v is not None:
                blf[ci[row[COL_CONC_STR]], bi[int(row[COL_BIN])]] = float(v)
    return CountTensors(clonotypes, concs[COL_CONC_STR].to_list(), concs[COL_CONC_VAL].to_numpy(), bins, R, cells, blf)


def _fit_chunk(
    args: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, tuple[float, float]],
) -> GridFit:
    return fit_titration(*args)


def _fit_parallel(
    R: np.ndarray,
    a: np.ndarray,
    conc: np.ndarray,
    log_values: np.ndarray,
    B: float,
    bounds: tuple[float, float],
) -> GridFit:
    V = R.shape[0]
    workers = min(available_cpus(), math.ceil(V / _CHUNK))
    if workers <= 1 or V < _PARALLEL_MIN_CLONOTYPES:
        return fit_titration(R, a, conc, log_values, B, bounds)
    tasks = [(R[s : s + _CHUNK], a[s : s + _CHUNK], conc, log_values, B, bounds) for s in range(0, V, _CHUNK)]
    log(f"  parallel fit: {workers} workers, {len(tasks)} chunks of <= {_CHUNK}")
    with ProcessPoolExecutor(max_workers=workers) as ex:
        parts = list(ex.map(_fit_chunk, tasks))
    return GridFit(
        np.concatenate([p.log10_kd for p in parts]),
        np.concatenate([p.amplitude for p in parts]),
        np.concatenate([p.log10_kd_sd for p in parts]),
        parts[0].log10_kd_grid,
    )


def _bin_log_fluorescence(t: CountTensors, params: FitParams) -> tuple[np.ndarray, bool]:
    """(C, B) log10 fluorescence of each bin at each concentration, and whether it was measured.

    Falls back to the fixed `bin_width * (label - 1)` scale, broadcast over concentrations, when no
    column is bound or the bound one cannot be used. Never refuses to run: a fixed scale is the
    behaviour every run had before the column existed.
    """
    fixed = np.broadcast_to(
        params.bin_width_decades * (t.bin_labels.astype(float) - 1.0),
        (t.concentrations.size, t.bin_labels.size),
    ).copy()
    m = t.bin_log_fluorescence
    if m is None:
        return fixed, False
    if not np.isfinite(m).all():
        log(
            "WARN: measured bin fluorescence is missing for some (concentration, bin) samples; "
            "falling back to the fixed bin scale for every sample."
        )
        return fixed, False
    if not (np.diff(m, axis=1) > 0).all():
        log(
            "WARN: measured bin fluorescence does not increase across bins at every concentration; "
            "the gates must be ordered low to high. Falling back to the fixed bin scale."
        )
        return fixed, False
    top_spread = float(m[:, -1].max() - m[:, -1].min())
    if top_spread > BIN_FLUORESCENCE_MAX_TOP_SPREAD:
        log(
            f"WARN: the measured top gate moves {top_spread:.2f} decades across concentrations "
            f"({m[:, -1].min():.2f} to {m[:, -1].max():.2f}). At the concentrations where it sits "
            f"low, a saturated binder's fluorescence is above the top gate, which the model cannot "
            f"represent, so the tightest binders would be unfittable and the ranking would invert. "
            f"Falling back to the fixed bin scale. Measured gate values need gates whose own mean "
            f"does not depend on which cells occupy them."
        )
        return fixed, False
    log(f"Bin fluorescence: measured per sample, range {m.min():.3f} to {m.max():.3f} (log10)")
    return m, True


def fit_poisson(
    reads: pl.DataFrame,
    signal_frame: pl.DataFrame,
    insufficient_map: dict[str, str | None],
    params: FitParams,
    cell_col: str | None,
    *,
    cells_absolute: bool,
    bin_fluorescence_col: str | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """cell_col: per-sample cell counts (cells_absolute=True) or sort fractions (False); None = uncorrected."""
    if not params.bin_width_decades > 0:
        raise InputValidationError("bin_width_decades must be > 0")
    if not 0 < params.min_plateau_fraction < 1:
        raise InputValidationError("min_plateau_fraction must be in (0, 1)")

    t = build_tensors(reads, cell_col, bin_fluorescence_col)
    a2, h = exposure(t.reads.sum(0), t.sample_cells, absolute=cells_absolute)
    R = t.reads * h[None]
    a = np.broadcast_to(a2, R.shape)
    log_values, measured_bins = _bin_log_fluorescence(t, params)

    zero = np.flatnonzero(t.concentrations == 0)
    if params.baseline_fluorescence is not None:
        B = float(params.baseline_fluorescence)
    elif zero.size:
        est = baseline(R[:, zero[0], :], a2[zero[0]], log_values[zero[0]], t.reads[:, zero[0], :])
        B = est.value
        log(
            f"0 M baseline computed from {est.n_kept} of {est.n_total} clonotypes "
            f"(>= {est.read_floor:.0f} raw 0 M reads)"
        )
        if est.contamination_ratio is not None and est.contamination_ratio > BASELINE_CONTAMINATION_RATIO:
            log(
                f"WARN: the 0 M clonotype set looks contaminated — the per-clonotype mean is "
                f"{est.contamination_ratio:.2f}x the pooled baseline, against 0.8-1.3 on a clean "
                f"library. Shallow sequencing-error clonotypes inflate the computed baseline, which "
                f"shifts every K_D. Check the clonotype set, or enter the 0 M baseline fluorescence."
            )
    else:
        B = None
    if B is None or not B > 0:
        raise InputValidationError(
            "poisson_log needs the 0 M baseline: add a no-antigen (0 M) sample or set baseline_fluorescence"
        )
    l_lo, l_hi = float(log_values.min()), float(log_values.max())
    top = 10.0**l_hi
    amp_lo = 10.0 ** (l_lo + params.min_plateau_fraction * (l_hi - l_lo)) - B
    if not amp_lo > 0:
        raise InputValidationError(
            f"min_plateau_fraction {params.min_plateau_fraction} puts the plateau below the 0 M baseline "
            f"(B = {B:.3g} bin-1 units); raise it"
        )
    log(
        f"poisson_log: B = {B:.4f}, bin width = {params.bin_width_decades} decades, "
        f"amplitude in [{amp_lo:.3g}, {top:.3g}]"
    )
    pos = t.concentrations[t.concentrations > 0]
    # Clonotypes the sufficiency gate rejected are reported Failed without a fit. Exposure and B above
    # still use the whole library.
    survivors = np.array([insufficient_map.get(key) is None for key in t.clonotypes]) & (pos.size > 0)
    n_fit = int(survivors.sum())
    log(f"Fitting {n_fit} clonotypes ({survivors.size - n_fit} pre-gated Failed)")
    fit = GridFit(*(np.full(survivors.size, np.nan) for _ in range(3)), np.empty(0))
    if n_fit:
        part = _fit_parallel(R[survivors], a[survivors], t.concentrations, log_values, B, (amp_lo, top))
        fit.log10_kd[survivors] = part.log10_kd
        fit.amplitude[survivors] = part.amplitude
        fit.log10_kd_sd[survivors] = part.log10_kd_sd
    lo, hi = (math.log10(pos.min()), math.log10(pos.max())) if pos.size else (math.nan, math.nan)
    conc_index = {c: i for i, c in enumerate(t.conc_strs)}
    signal = signal_frame.filter(pl.col(COL_CONC_VAL) > 0).sort([COL_CLONOTYPE, COL_CONC_VAL])
    by_clonotype = signal.partition_by(COL_CLONOTYPE, as_dict=True)

    rows: dict[str, list] = {k: [] for k in PER_CLONOTYPE_SCHEMA}
    fitted_parts = []
    for i, key in enumerate(t.clonotypes):
        reason = insufficient_map.get(key)
        if not survivors[i] and reason is None:
            reason = "insufficient_points"  # no non-zero concentration to fit
        lk = float(fit.log10_kd[i])
        censor = "below" if lk <= lo else "above" if lk >= hi else None
        if reason is not None:
            cls = "Failed"
        elif censor == "below":
            # Its own class: the clamped value is the floor, the smallest number in the table, so as
            # Good it would head any downstream ranking as an unorderable block of bounds.
            cls = "CensoredBelow"
        elif fit.log10_kd_sd[i] <= params.kd_sd_good:
            cls = "Good"
        else:
            cls = "Partial"
        # An above-range K_D is clamped and flagged, never failed. Failing it strips the weak-binder
        # tail out of the table and was the sole cause of every measured rank-agreement regression.
        has_fit = reason is None
        r2 = None
        sub = by_clonotype.get((key,))
        if has_fit and sub is not None and sub.height:
            sub_conc = sub[COL_CONC_VAL].to_numpy()
            if measured_bins:
                conc_rows = np.array([conc_index[c] for c in sub[COL_CONC_STR].to_list()])
                y_hat = predicted_mean_bin_measured(
                    lk, float(fit.amplitude[i]), sub_conc, B, log_values[conc_rows], t.bin_labels
                )
            else:
                y_hat = predicted_mean_bin(
                    lk,
                    float(fit.amplitude[i]),
                    sub_conc,
                    B,
                    params.bin_width_decades,
                    int(t.bin_labels[-1]),
                )
            ok = sub[SIGNAL].is_not_null().to_numpy()
            r2 = weighted_r2(sub[SIGNAL].to_numpy()[ok], y_hat[ok], sub[CLONOTYPE_READS_AT_CONC].to_numpy()[ok])
            r2 = None if r2 is None or math.isnan(r2) else r2
            if cls != "Failed":
                fitted_parts.append(
                    pl.DataFrame(
                        {COL_CLONOTYPE: [key] * sub.height, COL_CONC_STR: sub[COL_CONC_STR], "fittedMeanBin": y_hat},
                        schema=FITTED_MEAN_BIN_SCHEMA,
                    )
                )
        rows[COL_CLONOTYPE].append(key)
        rows["kd"].append(10.0 ** min(max(lk, lo), hi) if has_fit else None)
        rows["hillCoefficient"].append(1.0 if has_fit else None)
        rows["r2"].append(r2)
        rows["affinityClass"].append(cls)
        rows["fitFailureReason"].append(reason)
        rows["kdOutOfRange"].append(censor is not None if has_fit else None)
        rows["kdLog10Sd"].append(float(fit.log10_kd_sd[i]) if has_fit else None)
        rows["kdCensor"].append(censor if has_fit else None)
    per_clonotype = pl.DataFrame(rows, schema=PER_CLONOTYPE_SCHEMA)
    fitted = pl.concat(fitted_parts) if fitted_parts else pl.DataFrame(schema=FITTED_MEAN_BIN_SCHEMA)
    return per_clonotype, fitted
