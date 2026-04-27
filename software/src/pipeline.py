"""End-to-end fitting pipeline.

Orchestrates the stages that turn a reads frame into the three output frames:

  1. Load reads table (caller's responsibility).
  2. Validate R1-R5 and apply antigen filter (R4).
  3. Compute signal (R7 bin mode or R7b no-bin mode).
  4. Apply R8 weights + R9 floor filter.
  5. Compute global baseline B from c=0 points (R6).
  6. Detect hook effect per clonotype (R9b).
  7. Fit Hill equation per clonotype (R10) + weighted R² (R11).
  8. Classify (R12), flag kdOutOfRange (R14b).
  9. Assemble outputs: per-clonotype frame + meanBin + fittedMeanBin (R13, R14).

`run()` is the public entry point. The CLI in main.py wraps it with
argparse and file I/O.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np
import polars as pl

from classify import classify
from constants import (
    COL_ANTIGEN,
    COL_BIN,
    COL_CLONOTYPE,
    COL_CONC_STR,
    COL_CONC_VAL,
    DEFAULT_PARAMS,
    FitParams,
)
from hill_fit import fit_one_clonotype
from log import log
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
from normalization import SIGNAL, normalize
from output_build import (
    FITTED_MEAN_BIN_SCHEMA,
    PER_CLONOTYPE_SCHEMA,
    add_diagnostic_plot_columns,
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
            for done, fit in enumerate(
                executor.map(worker, xs, ys, ws, chunksize=chunksize), start=1
            ):
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
        kd_col[i] = fit.kd
        n_col[i] = fit.n
        r2_col[i] = fit.r2_w
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


def _build_outputs(
    per_clonotype: pl.DataFrame,
    fitted: pl.DataFrame,
    signal_frame: pl.DataFrame,
    reads: pl.DataFrame,
) -> dict[str, pl.DataFrame]:
    """Apply R14b kdOutOfRange flag and assemble the three-frame output dict."""
    min_max = (
        reads.filter(pl.col(COL_CONC_VAL) > 0)
        .select(
            pl.col(COL_CONC_VAL).min().alias("min"),
            pl.col(COL_CONC_VAL).max().alias("max"),
        )
        .row(0, named=True)
    )
    if min_max["min"] is not None:
        per_clonotype = flag_kd_out_of_range(per_clonotype, float(min_max["min"]), float(min_max["max"]))
    max_c = float(min_max["max"]) if min_max["max"] is not None else 1.0
    per_clonotype = add_diagnostic_plot_columns(per_clonotype, max_c)

    return {
        "per_clonotype": per_clonotype,
        "mean_bin": build_mean_bin_frame(signal_frame),
        "fitted_mean_bin": fitted,
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
