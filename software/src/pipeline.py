"""End-to-end fitting pipeline.

Orchestrates the stages that turn a reads frame into the three output frames. Validation, the
antigen filter and the signal are shared; `params.equation` then selects one of two fits, each with
its own gates, its own rating and its own out-of-range rule. Both emit the same columns.

Shared, steps 1-3:

  1. Load reads table (caller's responsibility).
  2. Validate R1-R5 and apply antigen filter (R4).
  3. Compute signal (R7 bin mode or R7b no-bin mode).

"hill_gate_index" — v2.4.0's fit, unchanged:

  4. Apply R8 weights + R9 floor filter, at the user's thresholds.
  5. Compute global baseline B from c=0 points (R6).
  6. Detect hook effect per clonotype (R9b).
  7. Fit the Hill equation per clonotype (R10) + weighted R² (R11).
  8. Classify on R² and the fitted exponent (R12), flag kdOutOfRange against the range (R14b).
  9. Assemble outputs: per-clonotype frame + meanBin + fittedMeanBin (R13, R14).

"log_fluorescence" — `_run_log_fluorescence`, bin mode only:

  4. Convert the mean gate index to mean log10 fluorescence; apply the floor and weights at the
     equation's own fixed thresholds.
  5. Fit log10(A*c/(c+K) + B) with the exponent fixed at 1 and the baseline free.
  6. Pin each K_D into the titrated range and record which wall, if any, set it.
  7. Classify on the bound and on SE(log10 K_D); kdOutOfRange reports the pin.
  8. Assemble the same outputs, with the fitted curve mapped back to the gate-index scale.

`run()` is the public entry point. The CLI in main.py wraps it with
argparse and file I/O.
"""

from __future__ import annotations

import math
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from functools import partial

import numpy as np
import polars as pl

from classify import classify, classify_log_fluorescence
from constants import (
    COL_ANTIGEN,
    COL_BIN,
    COL_CLONOTYPE,
    COL_CONC_STR,
    COL_CONC_VAL,
    DEFAULT_PARAMS,
    HILL_EXPONENT_FIXED,
    LF_MIN_CONCENTRATION_POINTS,
    LF_MIN_READS_PER_CONCENTRATION,
    FitParams,
)
from hill_fit import fit_one_clonotype
from io_layer import (
    apply_antigen_filter,
    canonicalize_concentration,
    max_bin_label,
    validate_antigen_filter,
    validate_bin_column,
    validate_bin_concentration_grid,
    validate_concentration_column,
    validate_reads_schema,
    validate_sample_metadata_uniqueness,
    validate_sort_fraction,
)
from log import log
from log_fluorescence import (
    GateScale,
    build_gate_scale,
    classify_curve,
    decades_to_gate_index,
    fit_one,
    gate_index_to_decades,
    pin_kd,
)
from normalization import SIGNAL, normalize
from output_build import (
    FITTED_MEAN_BIN_SCHEMA,
    PER_CLONOTYPE_SCHEMA,
    add_diagnostic_plot_columns,
    build_concentration_value_frame,
    build_mean_bin_frame,
    flag_kd_out_of_range,
)
from pre_fit import (
    WEIGHT,
    apply_floor_and_weights,
    classify_insufficient,
    compute_global_baseline,
    detect_hook_effect,
    split_c0,
)


def _validate_inputs(
    reads: pl.DataFrame,
    *,
    has_bin: bool,
    has_antigen: bool,
    target_antigen: str | None,
    antigen_column_ref: str | None,
) -> None:
    """Run all R1-R5 validators and emit accumulated warnings to stderr."""
    validate_reads_schema(reads, has_bin=has_bin, has_antigen=has_antigen)
    warnings: list[str] = []
    warnings += validate_concentration_column(reads, has_bin=has_bin)
    if has_bin:
        validate_bin_column(reads)
        warnings += validate_bin_concentration_grid(reads)
    warnings += validate_antigen_filter(reads, antigen_column_ref, target_antigen)
    validate_sample_metadata_uniqueness(reads, has_bin=has_bin, has_antigen=has_antigen)
    # Warnings go through the same log channel as progress messages. The Tengo workflow
    # captures the software's stdout via saveStdoutStream() and surfaces it in the block
    # Fit Log UI (main.tpl.tengo); stderr is not captured, so stderr writes would be invisible.
    for w in warnings:
        log(f"WARN: {w}")


# Pool start-up + pickling overhead dominates for small workloads; serial is faster
# until roughly this many clonotypes survive the pre-fit gates.
_PARALLEL_FIT_MIN_SURVIVORS = 50


def _execute_fits(
    worker,
    xs: list[np.ndarray],
    ys: list[np.ndarray],
    ws: list[np.ndarray],
    *,
    n_survivors: int,
    step: int,
) -> list:
    """Run worker(x, y, w) over each (x, y, w) triple; in-order results.

    Uses ProcessPoolExecutor when it pays off (worker count > 1 and enough
    survivors to amortize startup). The serial fallback preserves original
    semantics for tiny inputs and for environments where fork/spawn is a
    liability. Progress is logged from the parent as results arrive — `map`
    yields in input order, so we can count completions deterministically.
    """
    if n_survivors == 0:
        return []

    cpu_count = os.cpu_count() or 1
    max_workers = min(cpu_count, n_survivors)
    if max_workers > 1 and n_survivors >= _PARALLEL_FIT_MIN_SURVIVORS:
        # chunksize trades IPC overhead against load balancing. ~4 chunks per
        # worker keeps stragglers from stalling at the tail while batching
        # enough work per round-trip to hide pickle cost.
        chunksize = max(1, n_survivors // (max_workers * 4))
        log(f"  parallel fit: {max_workers} workers, chunksize={chunksize}")
        fits: list = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            for done, fit in enumerate(executor.map(worker, xs, ys, ws, chunksize=chunksize), start=1):
                fits.append(fit)
                if done % step == 0 or done == n_survivors:
                    log(f"  fitted {done}/{n_survivors} clonotypes")
        return fits

    fits = []
    for idx, (x, y, w) in enumerate(zip(xs, ys, ws)):
        fits.append(worker(x, y, w))
        done = idx + 1
        if done % step == 0 or done == n_survivors:
            log(f"  fitted {done}/{n_survivors} clonotypes")
    return fits


def _fit_all_clonotypes(
    fit_points: pl.DataFrame,
    *,
    all_clonotypes: list[str],
    insufficient_map: dict[str, str | None],
    hook_map: dict[str, bool],
    global_baseline: float | None,
    bin_mode: bool,
    max_bin: int | None,
    params: FitParams,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Run pre-fit gates, fit all surviving clonotypes, integrate results.

    Three phases:
      1. Gate — mark Failed/insufficient_*/non_monotonic_signal without fitting;
         collect (x, y, w) triples + fitted-row context for survivors.
      2. Execute — run fit_one_clonotype serially over the survivor list.
      3. Integrate — classify each FitResult and splat into pre-allocated columns.

    Partitions fit_points once (O(N)) instead of filtering per clonotype inside
    the loop. Output frames are built columnarly from python lists (nullable)
    and numpy chunks (fitted rows) — no per-clonotype dict churn.
    """
    sorted_points = fit_points.sort([COL_CLONOTYPE, COL_CONC_VAL])
    partitions = sorted_points.partition_by(COL_CLONOTYPE, as_dict=True)

    clonotypes_sorted = sorted(all_clonotypes)
    n_clon = len(clonotypes_sorted)

    # Python lists preserve None semantics end-to-end — numpy float arrays collapse
    # missing values into NaN, which polars treats as distinct from null (R13: null).
    kd_col: list[float | None] = [None] * n_clon
    n_col: list[float | None] = [None] * n_clon
    r2_col: list[float | None] = [None] * n_clon
    affinity_col: list[str] = [""] * n_clon
    reason_col: list[str | None] = [None] * n_clon

    # Phase 1 — gate clonotypes; build picklable task list for survivors.
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    ws: list[np.ndarray] = []
    # Per-task integration context — index, clonotype key, concentration strings, x.
    # Kept separate from task inputs so x is not shipped back over the pool.
    contexts: list[tuple[int, str, np.ndarray, np.ndarray]] = []

    for i, clonotype in enumerate(clonotypes_sorted):
        insufficient_reason = insufficient_map.get(clonotype)
        if insufficient_reason is not None:
            affinity_col[i] = "Failed"
            reason_col[i] = insufficient_reason
            continue
        if hook_map.get(clonotype, False):
            affinity_col[i] = "Failed"
            reason_col[i] = "non_monotonic_signal"
            continue

        sub = partitions[(clonotype,)]
        x = sub[COL_CONC_VAL].to_numpy()
        xs.append(x)
        ys.append(sub[SIGNAL].to_numpy())
        ws.append(sub[WEIGHT].to_numpy())
        contexts.append((i, clonotype, sub[COL_CONC_STR].to_numpy(), x))

    # Phase 2 — execute fits. partial binds the constant kwargs so the worker
    # receives one arg from each of (xs, ys, ws) per call.
    worker = partial(
        fit_one_clonotype,
        baseline_fixed=global_baseline,
        bin_mode=bin_mode,
        max_bin_label=max_bin,
    )
    n_survivors = len(xs)
    gated = n_clon - n_survivors
    log(f"Fitting {n_survivors} clonotypes ({gated} pre-gated Failed)")
    # Emit a progress line ~20 times over the fit loop (every 5%), with a floor of 1
    # so tiny inputs still report. Small overhead; stdout is flushed per line.
    step = max(1, n_survivors // 20)
    fits = _execute_fits(worker, xs, ys, ws, n_survivors=n_survivors, step=step)

    # Phase 3 — integrate results.
    fitted_keys: list[np.ndarray] = []
    fitted_conc_strs: list[np.ndarray] = []
    fitted_y_hats: list[np.ndarray] = []

    for fit, (i, clonotype, cs, x) in zip(fits, contexts):
        cls = classify(fit.r2_w, fit.n, fit.converged, params)
        # A signal with no variance leaves R² undefined. The curve the fitter returns is then
        # an artefact of the amplitude lower bound, not a measurement: on a flat bin mean it
        # reports a K_D of 28 M against a titration that stops at 1e-6. Report no numbers, as
        # for a fit that never converged.
        no_information = fit.r2_w is not None and not math.isfinite(fit.r2_w)
        kd_col[i] = None if no_information else fit.kd
        n_col[i] = None if no_information else fit.n
        r2_col[i] = None if no_information else fit.r2_w
        affinity_col[i] = cls.affinity_class
        reason_col[i] = cls.failure_reason

        if fit.converged and fit.y_hat is not None and cls.affinity_class != "Failed":
            k = x.shape[0]
            fitted_keys.append(np.full(k, clonotype, dtype=object))
            fitted_conc_strs.append(cs)
            fitted_y_hats.append(fit.y_hat)

    per_clonotype = pl.DataFrame(
        {
            COL_CLONOTYPE: clonotypes_sorted,
            "kd": kd_col,
            "hillCoefficient": n_col,
            "r2": r2_col,
            "affinityClass": affinity_col,
            "fitFailureReason": reason_col,
            "kdOutOfRange": [None] * n_clon,
        },
        schema=PER_CLONOTYPE_SCHEMA,
    )

    if fitted_keys:
        fitted = pl.DataFrame(
            {
                COL_CLONOTYPE: np.concatenate(fitted_keys),
                COL_CONC_STR: np.concatenate(fitted_conc_strs),
                "fittedMeanBin": np.concatenate(fitted_y_hats),
            },
            schema=FITTED_MEAN_BIN_SCHEMA,
        )
    else:
        fitted = pl.DataFrame({k: [] for k in FITTED_MEAN_BIN_SCHEMA}, schema=FITTED_MEAN_BIN_SCHEMA)

    return per_clonotype, fitted


def _concentration_bounds(reads: pl.DataFrame) -> tuple[float | None, float | None]:
    """Lowest and highest non-zero concentration, or (None, None) when there are none."""
    row = (
        reads.filter(pl.col(COL_CONC_VAL) > 0)
        .select(
            pl.col(COL_CONC_VAL).min().alias("min"),
            pl.col(COL_CONC_VAL).max().alias("max"),
        )
        .row(0, named=True)
    )
    lo = float(row["min"]) if row["min"] is not None else None
    hi = float(row["max"]) if row["max"] is not None else None
    return lo, hi


def _observed_statistics(floor_frame: pl.DataFrame, scale: GateScale) -> dict[str, tuple[float | None, float | None]]:
    """Per-clonotype observed rise and plateau height, both read by the bound test.

    Taken over every point that passed the read floor, the c=0 point included. `plateau_height` is
    the observed top as a fraction of the gate range — the statistic the level test compares to the
    plateau floor.
    """
    per_clonotype = floor_frame.group_by(COL_CLONOTYPE, maintain_order=True).agg(
        pl.col(SIGNAL).max().alias("top"),
        pl.col(SIGNAL).min().alias("bottom"),
    )
    stats = per_clonotype.with_columns(
        (pl.col("top") - pl.col("bottom")).alias("rise"),
        ((pl.col("top") - scale.floor) / scale.range_decades).alias("plateau_height"),
    )
    return {row[COL_CLONOTYPE]: (row["rise"], row["plateau_height"]) for row in stats.iter_rows(named=True)}


def _fit_all_log_fluorescence(
    fit_points: pl.DataFrame,
    *,
    all_clonotypes: list[str],
    insufficient_map: dict[str, str | None],
    stats_map: dict[str, tuple[float | None, float | None]],
    scale: GateScale,
    plateau_floor: float,
    log_lo: float,
    log_hi: float,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Fit every clonotype with "log_fluorescence", then pin, rate and integrate.

    Mirrors `_fit_all_clonotypes`'s three phases and its column-building, so both equations emit
    the same frames. `log_lo` and `log_hi` are log10 of the lowest and highest non-zero
    concentrations: the walls a K_D is pinned to.
    """
    sorted_points = fit_points.sort([COL_CLONOTYPE, COL_CONC_VAL])
    partitions = sorted_points.partition_by(COL_CLONOTYPE, as_dict=True)

    clonotypes_sorted = sorted(all_clonotypes)
    n_clon = len(clonotypes_sorted)

    kd_col: list[float | None] = [None] * n_clon
    n_col: list[float | None] = [None] * n_clon
    r2_col: list[float | None] = [None] * n_clon
    affinity_col: list[str] = [""] * n_clon
    reason_col: list[str | None] = [None] * n_clon
    out_of_range_col: list[bool | None] = [None] * n_clon

    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    ws: list[np.ndarray] = []
    contexts: list[tuple[int, str, np.ndarray, np.ndarray]] = []

    for i, clonotype in enumerate(clonotypes_sorted):
        insufficient_reason = insufficient_map.get(clonotype)
        if insufficient_reason is not None:
            affinity_col[i] = "Failed"
            reason_col[i] = insufficient_reason
            continue
        sub = partitions.get((clonotype,))
        if sub is None or sub.height == 0:
            affinity_col[i] = "Failed"
            reason_col[i] = "insufficient_reads"
            continue

        x = sub[COL_CONC_VAL].to_numpy()
        xs.append(x)
        ys.append(sub[SIGNAL].to_numpy())
        ws.append(sub[WEIGHT].to_numpy())
        contexts.append((i, clonotype, sub[COL_CONC_STR].to_numpy(), x))

    worker = partial(fit_one, scale=scale, plateau_floor=plateau_floor)
    n_survivors = len(xs)
    gated = n_clon - n_survivors
    log(f"Fitting {n_survivors} clonotypes ({gated} pre-gated Failed)")
    step = max(1, n_survivors // 20)
    fits = _execute_fits(worker, xs, ys, ws, n_survivors=n_survivors, step=step)

    fitted_keys: list[np.ndarray] = []
    fitted_conc_strs: list[np.ndarray] = []
    fitted_y_hats: list[np.ndarray] = []

    for fit, (i, clonotype, cs, x) in zip(fits, contexts):
        if not fit.converged or fit.log10_kd is None:
            affinity_col[i] = "Failed"
            reason_col[i] = "convergence_failure"
            continue

        rise, plateau_height = stats_map.get(clonotype, (None, None))
        curve_class = classify_curve(
            rise_decades=rise,
            plateau_height=plateau_height,
            plateau_floor=plateau_floor,
        )
        log10_kd, censor = pin_kd(fit.log10_kd, curve_class=curve_class, lo=log_lo, hi=log_hi)
        cls = classify_log_fluorescence(
            converged=True,
            censored=censor is not None,
            log10_kd_se=fit.log10_kd_se,
        )

        kd_col[i] = 10.0**log10_kd
        n_col[i] = HILL_EXPONENT_FIXED
        r2_col[i] = fit.r2
        affinity_col[i] = cls.affinity_class
        reason_col[i] = cls.failure_reason
        out_of_range_col[i] = censor is not None

        if fit.y_hat is not None and cls.affinity_class != "Failed":
            k = x.shape[0]
            fitted_keys.append(np.full(k, clonotype, dtype=object))
            fitted_conc_strs.append(cs)
            fitted_y_hats.append(decades_to_gate_index(fit.y_hat))

    per_clonotype = pl.DataFrame(
        {
            COL_CLONOTYPE: clonotypes_sorted,
            "kd": kd_col,
            "hillCoefficient": n_col,
            "r2": r2_col,
            "affinityClass": affinity_col,
            "fitFailureReason": reason_col,
            "kdOutOfRange": out_of_range_col,
        },
        schema=PER_CLONOTYPE_SCHEMA,
    )

    if fitted_keys:
        fitted = pl.DataFrame(
            {
                COL_CLONOTYPE: np.concatenate(fitted_keys),
                COL_CONC_STR: np.concatenate(fitted_conc_strs),
                "fittedMeanBin": np.concatenate(fitted_y_hats),
            },
            schema=FITTED_MEAN_BIN_SCHEMA,
        )
    else:
        fitted = pl.DataFrame({k: [] for k in FITTED_MEAN_BIN_SCHEMA}, schema=FITTED_MEAN_BIN_SCHEMA)

    return per_clonotype, fitted


def _run_log_fluorescence(
    signal_frame: pl.DataFrame,
    reads: pl.DataFrame,
    *,
    params: FitParams,
    has_bin: bool,
    max_bin: int | None,
) -> dict[str, pl.DataFrame]:
    """The "log_fluorescence" branch, from the mean-gate-index signal frame to the output dict."""
    if not has_bin:
        raise ValueError(
            "equation 'log_fluorescence' needs gates: the reads table has no bin column. "
            "No-bin mode is served by equation 'hill_gate_index'."
        )
    if not 0.0 < params.plateau_floor <= 1.0:
        raise ValueError(f"plateau_floor must be in (0, 1]; got {params.plateau_floor}")

    log_lo, log_hi = _concentration_bounds(reads)
    if log_lo is None or log_hi is None:
        raise ValueError("equation 'log_fluorescence' needs at least one non-zero concentration.")
    log_lo, log_hi = math.log10(log_lo), math.log10(log_hi)

    scale = build_gate_scale(max_bin)
    log(f"Gate scale: {scale.floor:.3f} to {scale.top:.3f} log10, {scale.range_decades:.3f} decades")
    log(f"Plateau floor: {params.plateau_floor:.4f} of the gate range")

    # The gate index and the log fluorescence differ by one affine map, so the observed signal is
    # converted rather than recomputed. Output frames stay on the gate-index scale.
    lf_frame = signal_frame.with_columns(gate_index_to_decades(pl.col(SIGNAL)).alias(SIGNAL))

    # The equation ignores the user's read floor and point count and uses its own constants.
    gate_params = replace(
        params,
        min_reads_per_concentration=LF_MIN_READS_PER_CONCENTRATION,
        min_concentration_points=LF_MIN_CONCENTRATION_POINTS,
    )
    log("Applying floor and weights")
    floor_frame = apply_floor_and_weights(lf_frame, gate_params)

    all_clonotypes = signal_frame[COL_CLONOTYPE].unique().to_list()
    log(f"Classifying {len(all_clonotypes)} clonotypes for sufficiency")
    insufficient = classify_insufficient(floor_frame, all_clonotypes, gate_params)
    insufficient_map = dict(zip(insufficient[COL_CLONOTYPE].to_list(), insufficient["insufficient_reason"].to_list()))

    stats_map = _observed_statistics(floor_frame, scale)
    _, fit_points = split_c0(floor_frame)

    per_clonotype, fitted = _fit_all_log_fluorescence(
        fit_points,
        all_clonotypes=all_clonotypes,
        insufficient_map=insufficient_map,
        stats_map=stats_map,
        scale=scale,
        plateau_floor=params.plateau_floor,
        log_lo=log_lo,
        log_hi=log_hi,
    )

    log("Assembling output frames")
    outputs = _build_outputs(per_clonotype, fitted, signal_frame, reads, flag_out_of_range=False)
    log(
        f"Pipeline done: per_clonotype={outputs['per_clonotype'].height}, "
        f"mean_bin={outputs['mean_bin'].height}, "
        f"fitted_mean_bin={outputs['fitted_mean_bin'].height}"
    )
    return outputs


def _build_outputs(
    per_clonotype: pl.DataFrame,
    fitted: pl.DataFrame,
    signal_frame: pl.DataFrame,
    reads: pl.DataFrame,
    *,
    flag_out_of_range: bool = True,
) -> dict[str, pl.DataFrame]:
    """Apply R14b kdOutOfRange flag and assemble the four-frame output dict.

    `flag_out_of_range=False` keeps the flag the caller already set. "log_fluorescence" clamps its
    K_D into the titrated range, so a comparison against that range would read every pinned value
    as in-range.
    """
    min_c, max_c_opt = _concentration_bounds(reads)
    if flag_out_of_range and min_c is not None and max_c_opt is not None:
        per_clonotype = flag_kd_out_of_range(per_clonotype, min_c, max_c_opt)
    max_c = max_c_opt if max_c_opt is not None else 1.0
    per_clonotype = add_diagnostic_plot_columns(per_clonotype, max_c)

    # Don't add a .sort() here. Phase 3 builds both frames in deterministic order
    # (clonotypes_sorted; contexts iterated in conc_val-numeric order from
    # sorted_points). Sorting again rewrites the bytes and breaks dedup against
    # cached resources from earlier runs.
    return {
        "per_clonotype": per_clonotype,
        "mean_bin": build_mean_bin_frame(signal_frame),
        "fitted_mean_bin": fitted,
        "concentration_value": build_concentration_value_frame(signal_frame),
    }


def run(
    reads: pl.DataFrame,
    *,
    params: FitParams = DEFAULT_PARAMS,
    target_antigen: str | None = None,
    antigen_column_ref: str | None = None,
    sort_fraction_col: str | None = None,
) -> dict[str, pl.DataFrame]:
    """Execute the full pipeline on an in-memory reads frame.

    `sort_fraction_col`, when set, activates the Adams, Mora, Walczak, Kinney 2016
    eq. A3 correction inside `compute_mean_bin`. The column must be present in
    `reads`; validation runs before normalization so a bad upstream join is caught
    before any downstream computation. Only consulted in bin mode.
    """
    log(f"Pipeline start: {reads.height} reads rows")
    reads = canonicalize_concentration(reads)
    has_bin = COL_BIN in reads.columns
    has_antigen = COL_ANTIGEN in reads.columns
    log(f"Mode: bin={has_bin}, antigen_column={has_antigen}")

    log("Validating inputs (R1-R5)")
    _validate_inputs(
        reads,
        has_bin=has_bin,
        has_antigen=has_antigen,
        target_antigen=target_antigen,
        antigen_column_ref=antigen_column_ref,
    )

    if sort_fraction_col is not None and has_bin:
        validate_sort_fraction(reads, sort_fraction_col)
        n_conc = reads.select(COL_CONC_STR).n_unique()
        log(f"sort_fraction validated: {n_conc} concentrations, all sums within 1e-3 of 1.0")
        log(f"Mean-bin correction: FACS-weighted (column '{sort_fraction_col}')")
    elif has_bin:
        log("Mean-bin correction: uncorrected")

    if target_antigen is not None:
        log(f"Applying antigen filter: {target_antigen}")
    reads = apply_antigen_filter(reads, target_antigen)
    mbl = max_bin_label(reads) if has_bin else None

    log("Normalizing signals")
    signal_frame = normalize(
        reads,
        bin_mode=has_bin,
        sort_fraction_col=sort_fraction_col if has_bin else None,
    )

    log(f"Equation: {params.equation}")
    if params.equation == "log_fluorescence":
        return _run_log_fluorescence(signal_frame, reads, params=params, has_bin=has_bin, max_bin=mbl)

    log("Applying floor and weights")
    floor_frame = apply_floor_and_weights(signal_frame, params)

    all_clonotypes = signal_frame[COL_CLONOTYPE].unique().to_list()
    log(f"Classifying {len(all_clonotypes)} clonotypes for sufficiency")
    insufficient = classify_insufficient(floor_frame, all_clonotypes, params)
    insufficient_map = dict(zip(insufficient[COL_CLONOTYPE].to_list(), insufficient["insufficient_reason"].to_list()))

    c0_points, fit_points = split_c0(floor_frame)
    global_b = compute_global_baseline(c0_points)
    b_str = f"{global_b:.4f}" if global_b is not None else "none"
    log(f"Global baseline B = {b_str} (c=0 points: {c0_points.height})")

    log("Detecting hook effect")
    hook = detect_hook_effect(fit_points, bin_mode=has_bin, params=params)
    hook_map = dict(zip(hook[COL_CLONOTYPE].to_list(), hook["hook_flag"].to_list()))

    per_clonotype, fitted = _fit_all_clonotypes(
        fit_points,
        all_clonotypes=all_clonotypes,
        insufficient_map=insufficient_map,
        hook_map=hook_map,
        global_baseline=global_b,
        bin_mode=has_bin,
        max_bin=mbl,
        params=params,
    )

    log("Assembling output frames")
    outputs = _build_outputs(per_clonotype, fitted, signal_frame, reads)
    log(
        f"Pipeline done: per_clonotype={outputs['per_clonotype'].height}, "
        f"mean_bin={outputs['mean_bin'].height}, "
        f"fitted_mean_bin={outputs['fitted_mean_bin'].height}"
    )
    return outputs
