"""Bin-mode pipeline for the poisson_log model (Adams, Mora, Walczak, Kinney 2016, Appendix E).

Turns the long reads frame into count tensors, fits every clonotype with `poisson_fit`, and returns
the per-clonotype and fitted-mean-bin frames in the same schemas as the hill_free_n path.
"""

from __future__ import annotations

import math
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
import polars as pl

from constants import COL_BIN, COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL, COL_READS, COL_SAMPLE, FitParams
from hill_fit import weighted_r2
from io_layer import InputValidationError
from log import log
from normalization import CLONOTYPE_READS_AT_CONC, SIGNAL
from output_build import FITTED_MEAN_BIN_SCHEMA, PER_CLONOTYPE_SCHEMA
from poisson_fit import GridFit, baseline, exposure, fit_titration, predicted_mean_bin

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


def build_tensors(reads: pl.DataFrame, cell_count_col: str | None) -> CountTensors:
    clonotypes = sorted(reads[COL_CLONOTYPE].unique().to_list())
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
    return CountTensors(clonotypes, concs[COL_CONC_STR].to_list(), concs[COL_CONC_VAL].to_numpy(), bins, R, cells)


def _fit_chunk(args) -> GridFit:
    return fit_titration(*args)


def _fit_parallel(R, a, conc, log_values, B, bounds) -> GridFit:
    V = R.shape[0]
    workers = min(os.cpu_count() or 1, math.ceil(V / _CHUNK))
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


def fit_poisson(
    reads: pl.DataFrame,
    signal_frame: pl.DataFrame,
    insufficient_map: dict[str, str | None],
    params: FitParams,
    cell_col: str | None,
    *,
    cells_absolute: bool,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """cell_col: per-sample cell counts (cells_absolute=True) or sort fractions (False); None = uncorrected."""
    if not params.bin_width_decades > 0:
        raise InputValidationError("bin_width_decades must be > 0")
    if not 0 < params.min_plateau_fraction < 1:
        raise InputValidationError("min_plateau_fraction must be in (0, 1)")

    t = build_tensors(reads, cell_col)
    a2, h = exposure(t.reads.sum(0), t.sample_cells, absolute=cells_absolute)
    R = t.reads * h[None]
    a = np.broadcast_to(a2, R.shape)
    log_values = params.bin_width_decades * (t.bin_labels.astype(float) - 1.0)

    zero = np.flatnonzero(t.concentrations == 0)
    if params.baseline_fluorescence is not None:
        B = float(params.baseline_fluorescence)
    elif zero.size:
        B = baseline(R[:, zero[0], :], a2[zero[0]], log_values)
    else:
        B = None
    if B is None or not B > 0:
        raise InputValidationError(
            "poisson_log needs the 0 M baseline: add a no-antigen (0 M) sample or set baseline_fluorescence"
        )
    top = 10.0 ** log_values[-1]
    amp_lo = 10.0 ** (params.min_plateau_fraction * log_values[-1]) - B
    if not amp_lo > 0:
        raise InputValidationError(
            f"min_plateau_fraction {params.min_plateau_fraction} puts the plateau below the 0 M baseline "
            f"(B = {B:.3g} bin-1 units); raise it"
        )
    log(
        f"poisson_log: B = {B:.4f}, bin width = {params.bin_width_decades} decades, "
        f"amplitude in [{amp_lo:.3g}, {top:.3g}]"
    )
    log(f"Fitting {R.shape[0]} clonotypes")
    fit = _fit_parallel(R, a, t.concentrations, log_values, B, (amp_lo, top))

    pos = t.concentrations[t.concentrations > 0]
    lo, hi = math.log10(pos.min()), math.log10(pos.max())
    signal = signal_frame.filter(pl.col(COL_CONC_VAL) > 0).sort([COL_CLONOTYPE, COL_CONC_VAL])
    by_clonotype = signal.partition_by(COL_CLONOTYPE, as_dict=True)

    rows: dict[str, list] = {k: [] for k in PER_CLONOTYPE_SCHEMA}
    fitted_parts = []
    for i, key in enumerate(t.clonotypes):
        reason = insufficient_map.get(key)
        lk = float(fit.log10_kd[i])
        censor = "below" if lk <= lo else "above" if lk >= hi else None
        sub = by_clonotype.get((key,))
        y_hat = None
        if reason is None and sub is not None and sub.height:
            y_hat = predicted_mean_bin(
                lk,
                float(fit.amplitude[i]),
                sub[COL_CONC_VAL].to_numpy(),
                B,
                params.bin_width_decades,
                int(t.bin_labels[-1]),
            )
        if reason is not None:
            cls = "Failed"
        elif censor == "above":
            cls, reason = "Failed", "kd_above_range"
        elif censor == "below" or fit.log10_kd_sd[i] <= params.kd_sd_good:
            cls = "Good"
        else:
            cls = "Partial"
        fitted = reason is None or reason == "kd_above_range"
        rows[COL_CLONOTYPE].append(key)
        rows["kd"].append(10.0 ** min(max(lk, lo), hi) if fitted else None)
        rows["hillCoefficient"].append(1.0 if fitted else None)
        r2 = None
        if y_hat is not None:
            ok = sub[SIGNAL].is_not_null().to_numpy()
            r2 = weighted_r2(sub[SIGNAL].to_numpy()[ok], y_hat[ok], sub[CLONOTYPE_READS_AT_CONC].to_numpy()[ok])
            r2 = None if r2 is None or math.isnan(r2) else r2
        rows["r2"].append(r2)
        rows["affinityClass"].append(cls)
        rows["fitFailureReason"].append(reason)
        rows["kdOutOfRange"].append(censor is not None if fitted else None)
        rows["kdLog10Sd"].append(float(fit.log10_kd_sd[i]) if fitted else None)
        rows["kdCensor"].append(censor if fitted else None)
        if cls != "Failed" and y_hat is not None:
            fitted_parts.append(
                pl.DataFrame(
                    {COL_CLONOTYPE: [key] * sub.height, COL_CONC_STR: sub[COL_CONC_STR], "fittedMeanBin": y_hat},
                    schema=FITTED_MEAN_BIN_SCHEMA,
                )
            )
    per_clonotype = pl.DataFrame(rows, schema=PER_CLONOTYPE_SCHEMA)
    fitted = (
        pl.concat(fitted_parts)
        if fitted_parts
        else pl.DataFrame({k: [] for k in FITTED_MEAN_BIN_SCHEMA}, schema=FITTED_MEAN_BIN_SCHEMA)
    )
    return per_clonotype, fitted
