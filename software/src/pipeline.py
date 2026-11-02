"""End-to-end fitting pipeline.

Orchestrates the stages that turn a reads frame into the three output frames. Validation, the
antigen filter and the signal are shared; `params.equation` then selects one of two fits, each with
its own gates, its own rating and its own out-of-range rule. Both emit the same columns.

Shared, steps 1-3:

  1. Load reads table (caller's responsibility).
  2. Validate R1-R5 and apply antigen filter (R4).
  3. Compute signal (R7 bin mode or R7b no-bin mode).

"hill_gate_index" — v2.4.0's fit, plus the undefined-R² guard (see classify.r2_undefined):

  4. Apply R8 weights + R9 floor filter, at the user's thresholds.
  5. Compute global baseline B from c=0 points (R6).
  6. Detect hook effect per clonotype (R9b).
  7. Fit the Hill equation per clonotype (R10) + weighted R² (R11).
  8. Classify on R² and the fitted exponent (R12), flag kdOutOfRange against the range (R14b).
  9. Assemble outputs: per-clonotype frame + meanBin + fittedMeanBin (R13, R14).

"log_fluorescence" — `_run_log_fluorescence`, bin mode only:

  4. Convert the mean gate index to mean log10 fluorescence; apply the floor and weights at the
     equation's own fixed thresholds.
  5. Fit log10(A*c/(c+K) + B) with the exponent fixed at 1. With 0 M samples every clonotype
     shares one baseline, the read-weighted mean 0 M signal; without them each fits its own.
  6. Pin each K_D into the titrated range and record which wall, if any, set it.
  7. Classify on the bound and on SE(log10 K_D); kdOutOfRange reports the pin.
  8. Assemble the same outputs, with the fitted curve on the observed signal's scale: mean log10
     fluorescence with measured gates, the mean gate index with assumed ones.

`run()` is the public entry point. The CLI in main.py wraps it with
argparse and file I/O.
"""

from __future__ import annotations

import math
import os
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from functools import partial
from pathlib import Path

import numpy as np
import polars as pl

from classify import Classification, classify, classify_log_fluorescence, r2_undefined
from constants import (
    COL_ANTIGEN,
    COL_BIN,
    COL_BIN_FLUORESCENCE,
    COL_CLONOTYPE,
    COL_CONC_STR,
    COL_CONC_VAL,
    COL_READS,
    DEFAULT_PARAMS,
    GATE_WIDTH_DECADES,
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
    validate_bin_fluorescence,
    validate_concentration_column,
    validate_reads_schema,
    validate_sample_metadata_uniqueness,
    validate_sort_fraction,
)
from log import log
from log_fluorescence import (
    GateScale,
    build_gate_scale,
    build_measured_gate_scale,
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
    compute_read_weighted_baseline,
    detect_hook_effect,
    split_c0,
)
from synonymous import SynonymousInputs, combine_synonymous


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


def _host_cpus() -> int:
    """CPUs this process may be scheduled on: the affinity set where the OS reports one."""
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def _cgroup_cpu_quota(cgroup_root: Path) -> int | None:
    """The container's CPU limit in whole CPUs, rounded up; None when unlimited or unreadable.

    Reads the cgroup mounted at `cgroup_root` only. A malformed or unreadable file falls back to
    the host count rather than failing the fit.
    """
    try:
        v2 = cgroup_root / "cpu.max"
        if v2.is_file():
            quota, period = v2.read_text().split()
            if quota == "max":
                return None
            return math.ceil(int(quota) / int(period))
        v1_quota = cgroup_root / "cpu" / "cpu.cfs_quota_us"
        v1_period = cgroup_root / "cpu" / "cpu.cfs_period_us"
        if v1_quota.is_file() and v1_period.is_file():
            quota = int(v1_quota.read_text())
            if quota <= 0:
                return None
            return math.ceil(quota / int(v1_period.read_text()))
    except (OSError, ValueError, ZeroDivisionError):
        return None
    return None


def available_cpus(cgroup_root: Path = Path("/sys/fs/cgroup")) -> int:
    """CPUs the job was granted. In a container os.cpu_count() reports the host's cores."""
    host = _host_cpus()
    quota = _cgroup_cpu_quota(cgroup_root)
    return max(1, host if quota is None else min(host, quota))


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

    max_workers = min(available_cpus(), n_survivors)
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


# Per-clonotype output columns after the key, in PER_CLONOTYPE_SCHEMA order.
_PER_CLONOTYPE_FIELDS = (
    "kd",
    "hillCoefficient",
    "r2",
    "affinityClass",
    "fitFailureReason",
    "kdOutOfRange",
    "kdLog10Se",
)


class _FitBatch:
    """One equation's pass over every clonotype: gate, fit the survivors, collect both frames.

    Both equations run through this. Each supplies its gate reasons, its worker and its own
    integration of a fit result, so they emit the same columns in the same row order.

    Gating partitions `fit_points` once (O(N)) instead of filtering per clonotype. Columns are
    Python lists so None stays null end to end: a numpy float array turns a missing value into
    NaN, which polars treats as distinct from null (R13).
    """

    def __init__(
        self,
        fit_points: pl.DataFrame,
        all_clonotypes: list[str],
        gate_reason: Callable[[str], str | None],
    ) -> None:
        partitions = fit_points.sort([COL_CLONOTYPE, COL_CONC_VAL]).partition_by(COL_CLONOTYPE, as_dict=True)
        self.clonotypes = sorted(all_clonotypes)
        n_clon = len(self.clonotypes)
        self.columns: dict[str, list] = {name: [None] * n_clon for name in _PER_CLONOTYPE_FIELDS}
        self.columns["affinityClass"] = [""] * n_clon

        self._xs: list[np.ndarray] = []
        self._ys: list[np.ndarray] = []
        self._ws: list[np.ndarray] = []
        # Per-task integration context: output index, clonotype key, concentration strings. Kept
        # apart from the task inputs so nothing beyond (x, y, w) is shipped over the pool.
        self._contexts: list[tuple[int, str, np.ndarray]] = []
        self._fitted_keys: list[np.ndarray] = []
        self._fitted_conc_strs: list[np.ndarray] = []
        self._fitted_curves: list[np.ndarray] = []

        for i, clonotype in enumerate(self.clonotypes):
            reason = gate_reason(clonotype)
            if reason is not None:
                self.fail(i, reason)
                continue
            # A clonotype that passed the gates has at least one fit point, so its partition exists.
            sub = partitions[(clonotype,)]
            self._xs.append(sub[COL_CONC_VAL].to_numpy())
            self._ys.append(sub[SIGNAL].to_numpy())
            self._ws.append(sub[WEIGHT].to_numpy())
            self._contexts.append((i, clonotype, sub[COL_CONC_STR].to_numpy()))

    def fit(self, worker) -> list[tuple[int, str, np.ndarray, object]]:
        """Run `worker` over the survivors; returns (index, clonotype, concentration strings, fit)."""
        n_survivors = len(self._xs)
        log(f"Fitting {n_survivors} clonotypes ({len(self.clonotypes) - n_survivors} pre-gated Failed)")
        # A progress line about every 5 percent of the fits, and at least every fit for tiny inputs.
        step = max(1, n_survivors // 20)
        fits = _execute_fits(worker, self._xs, self._ys, self._ws, n_survivors=n_survivors, step=step)
        return [(i, clonotype, cs, fit) for (i, clonotype, cs), fit in zip(self._contexts, fits)]

    def fail(self, i: int, reason: str) -> None:
        self.columns["affinityClass"][i] = "Failed"
        self.columns["fitFailureReason"][i] = reason

    def rate(self, i: int, cls: Classification) -> None:
        self.columns["affinityClass"][i] = cls.affinity_class
        self.columns["fitFailureReason"][i] = cls.failure_reason

    def add_curve(self, clonotype: str, conc_strs: np.ndarray, curve: np.ndarray) -> None:
        self._fitted_keys.append(np.full(conc_strs.shape[0], clonotype, dtype=object))
        self._fitted_conc_strs.append(conc_strs)
        self._fitted_curves.append(curve)

    def frames(self) -> tuple[pl.DataFrame, pl.DataFrame]:
        per_clonotype = pl.DataFrame({COL_CLONOTYPE: self.clonotypes, **self.columns}, schema=PER_CLONOTYPE_SCHEMA)
        if not self._fitted_keys:
            return per_clonotype, pl.DataFrame({k: [] for k in FITTED_MEAN_BIN_SCHEMA}, schema=FITTED_MEAN_BIN_SCHEMA)
        fitted = pl.DataFrame(
            {
                COL_CLONOTYPE: np.concatenate(self._fitted_keys),
                COL_CONC_STR: np.concatenate(self._fitted_conc_strs),
                "fittedMeanBin": np.concatenate(self._fitted_curves),
            },
            schema=FITTED_MEAN_BIN_SCHEMA,
        )
        return per_clonotype, fitted


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
    """ "hill_gate_index": gate, fit and classify every clonotype.

    The gates mark Failed without fitting: insufficient reads or points first, then the hook effect.
    """

    def gate_reason(clonotype: str) -> str | None:
        reason = insufficient_map.get(clonotype)
        if reason is None and hook_map.get(clonotype, False):
            return "non_monotonic_signal"
        return reason

    batch = _FitBatch(fit_points, all_clonotypes, gate_reason)
    worker = partial(fit_one_clonotype, baseline_fixed=global_baseline, bin_mode=bin_mode, max_bin_label=max_bin)
    columns = batch.columns
    for i, clonotype, conc_strs, fit in batch.fit(worker):
        cls = classify(fit.r2_w, fit.n, fit.converged, params)
        # An undefined R² rates Failed, and its fitted numbers are an artefact; see r2_undefined.
        if not r2_undefined(fit.r2_w):
            columns["kd"][i] = fit.kd
            columns["hillCoefficient"][i] = fit.n
            columns["r2"][i] = fit.r2_w
        batch.rate(i, cls)
        if cls.affinity_class != "Failed":
            batch.add_curve(clonotype, conc_strs, fit.y_hat)
    return batch.frames()


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
    log10_conc_lo: float,
    log10_conc_hi: float,
    baseline: float | None,
    gates_measured: bool,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """ "log_fluorescence": fit every clonotype, then pin, rate and integrate.

    `log10_conc_lo` and `log10_conc_hi` are log10 of the lowest and highest non-zero
    concentrations: the walls a K_D is pinned to.

    The fitted curve is written on the observed signal's scale. With measured gates both are mean
    log10 fluorescence. With assumed gates the observed signal is the mean gate index, so the curve
    is mapped back to it.
    """
    batch = _FitBatch(fit_points, all_clonotypes, insufficient_map.get)
    worker = partial(fit_one, scale=scale, plateau_floor=plateau_floor, baseline=baseline)
    columns = batch.columns
    for i, clonotype, conc_strs, fit in batch.fit(worker):
        if not fit.converged:
            batch.fail(i, "convergence_failure")
            continue

        rise, plateau_height = stats_map[clonotype]
        curve_class = classify_curve(rise_decades=rise, plateau_height=plateau_height, plateau_floor=plateau_floor)
        log10_kd, censor = pin_kd(fit.log10_kd, curve_class=curve_class, lo=log10_conc_lo, hi=log10_conc_hi)

        columns["kd"][i] = 10.0**log10_kd
        columns["hillCoefficient"][i] = HILL_EXPONENT_FIXED
        columns["r2"][i] = fit.r2
        columns["kdOutOfRange"][i] = censor is not None
        # A pinned K_D is a bound. The fit's SE describes the unpinned estimate, so a bound gets none.
        se = fit.log10_kd_se
        columns["kdLog10Se"][i] = se if censor is None and se is not None and math.isfinite(se) else None
        batch.rate(i, classify_log_fluorescence(censored=censor is not None, log10_kd_se=se))
        batch.add_curve(clonotype, conc_strs, fit.y_hat if gates_measured else decades_to_gate_index(fit.y_hat))
    return batch.frames()


def _run_hill_gate_index(
    signal_frame: pl.DataFrame,
    *,
    params: FitParams,
    has_bin: bool,
    max_bin: int | None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The "hill_gate_index" branch, from the signal frame to the per-clonotype and fitted frames."""
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

    return _fit_all_clonotypes(
        fit_points,
        all_clonotypes=all_clonotypes,
        insufficient_map=insufficient_map,
        hook_map=hook_map,
        global_baseline=global_b,
        bin_mode=has_bin,
        max_bin=max_bin,
        params=params,
    )


def _run_log_fluorescence(
    signal_frame: pl.DataFrame,
    reads: pl.DataFrame,
    *,
    params: FitParams,
    has_bin: bool,
    max_bin: int | None,
    conc_bounds: tuple[float | None, float | None],
    gates_measured: bool,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The "log_fluorescence" branch, from the signal frame to the per-clonotype and fitted frames.

    `gates_measured` says the signal frame already carries mean log10 fluorescence, because the
    caller passed measured gate positions to `normalize`. The affine conversion below is then
    skipped, and the gate scale comes from the measured positions rather than the gate labels.
    """
    if not has_bin:
        raise ValueError(
            "equation 'log_fluorescence' needs gates: the reads table has no bin column. "
            "No-bin mode is served by equation 'hill_gate_index'."
        )
    if not 0.0 < params.plateau_floor <= 1.0:
        raise ValueError(f"plateau_floor must be in (0, 1]; got {params.plateau_floor}")

    conc_lo, conc_hi = conc_bounds
    if conc_lo is None or conc_hi is None:
        raise ValueError("equation 'log_fluorescence' needs at least one non-zero concentration.")

    if gates_measured:
        scale = build_measured_gate_scale(reads)
        lf_frame = signal_frame
    else:
        scale = build_gate_scale(max_bin)
        # The gate index and the log fluorescence differ by one affine map, so the observed
        # signal is converted rather than recomputed.
        lf_frame = signal_frame.with_columns(gate_index_to_decades(pl.col(SIGNAL)).alias(SIGNAL))
    log(f"Gate scale: {scale.floor:.3f} to {scale.top:.3f} log10, {scale.range_decades:.3f} decades")
    log(f"Plateau floor: {params.plateau_floor:.4f} of the gate range")

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
    # One baseline for every clonotype: the read-weighted mean 0 M signal. Without 0 M samples,
    # each clonotype fits its own.
    c0_points, fit_points = split_c0(floor_frame)
    baseline = compute_read_weighted_baseline(c0_points)
    if baseline is None:
        log("Baseline: fitted per clonotype (no 0 M samples)")
    else:
        log(f"Baseline: {baseline:.4f} log10, the read-weighted mean 0 M signal of {c0_points.height} clonotypes")

    return _fit_all_log_fluorescence(
        fit_points,
        all_clonotypes=all_clonotypes,
        insufficient_map=insufficient_map,
        stats_map=stats_map,
        scale=scale,
        plateau_floor=params.plateau_floor,
        log10_conc_lo=math.log10(conc_lo),
        log10_conc_hi=math.log10(conc_hi),
        baseline=baseline,
        gates_measured=gates_measured,
    )


def _build_outputs(
    per_clonotype: pl.DataFrame,
    fitted: pl.DataFrame,
    signal_frame: pl.DataFrame,
    conc_bounds: tuple[float | None, float | None],
    *,
    flag_out_of_range: bool,
) -> dict[str, pl.DataFrame]:
    """Apply R14b kdOutOfRange flag and assemble the four-frame output dict.

    `flag_out_of_range=False` keeps the flag the caller already set. "log_fluorescence" clamps its
    K_D into the titrated range, so a comparison against that range would read every pinned value
    as in-range.
    """
    min_c, max_c = conc_bounds
    if flag_out_of_range and min_c is not None and max_c is not None:
        per_clonotype = flag_kd_out_of_range(per_clonotype, min_c, max_c)
    per_clonotype = add_diagnostic_plot_columns(per_clonotype, max_c if max_c is not None else 1.0)

    # Don't add a .sort() here. _FitBatch builds both frames in deterministic order (sorted
    # clonotypes; fit points in numeric concentration order). Sorting again rewrites the bytes and
    # breaks dedup against cached resources from earlier runs.
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
    synonymous: SynonymousInputs | None = None,
) -> dict:
    """Execute the full pipeline on an in-memory reads frame.

    `sort_fraction_col`, when set, activates the Adams, Mora, Walczak, Kinney 2016
    eq. A3 correction inside `compute_mean_bin`. The column must be present in
    `reads`; validation runs before normalization so a bad upstream join is caught
    before any downstream computation. Only consulted in bin mode.

    `synonymous`, when set, adds `per_variant` (one row per amino-acid variant) and
    `synonymous_summary` (a JSON-ready dict) to the outputs. Every other output is unchanged.
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

    if target_antigen is not None:
        log(f"Applying antigen filter: {target_antigen}")
    reads = apply_antigen_filter(reads, target_antigen)

    # Sort fractions sum to 1 within one antigen's sorts, so the check reads the filtered frame.
    if sort_fraction_col is not None and has_bin:
        validate_sort_fraction(reads, sort_fraction_col)
        n_conc = reads.select(COL_CONC_STR).n_unique()
        log(f"sort_fraction validated: {n_conc} concentrations, all sums within 1e-3 of 1.0")
        log(f"Mean-bin correction: FACS-weighted (column '{sort_fraction_col}')")
    elif has_bin:
        log("Mean-bin correction: uncorrected")
    mbl = max_bin_label(reads) if has_bin else None

    # Measured gate positions are read by "log_fluorescence" only, and only in bin mode. They must
    # enter here, before the gates collapse into one number a concentration, because a position
    # varies per gate per concentration.
    measured_gates = (
        COL_BIN_FLUORESCENCE
        if has_bin and params.equation == "log_fluorescence" and COL_BIN_FLUORESCENCE in reads.columns
        else None
    )
    if measured_gates is not None:
        validate_bin_fluorescence(reads)
        log(f"Gate positions: measured (column '{COL_BIN_FLUORESCENCE}')")
    elif params.equation == "log_fluorescence":
        log(f"Gate positions: assumed, {GATE_WIDTH_DECADES} decades a gate")

    log("Normalizing signals")
    signal_frame = normalize(
        reads,
        bin_mode=has_bin,
        sort_fraction_col=sort_fraction_col if has_bin else None,
        gate_position_col=measured_gates,
    )

    conc_bounds = _concentration_bounds(reads)
    log(f"Equation: {params.equation}")
    if params.equation == "log_fluorescence":
        per_clonotype, fitted = _run_log_fluorescence(
            signal_frame,
            reads,
            params=params,
            has_bin=has_bin,
            max_bin=mbl,
            conc_bounds=conc_bounds,
            gates_measured=measured_gates is not None,
        )
    else:
        per_clonotype, fitted = _run_hill_gate_index(signal_frame, params=params, has_bin=has_bin, max_bin=mbl)

    log("Assembling output frames")
    outputs = _build_outputs(
        per_clonotype,
        fitted,
        signal_frame,
        conc_bounds,
        flag_out_of_range=params.equation == "hill_gate_index",
    )
    log(
        f"Pipeline done: per_clonotype={outputs['per_clonotype'].height}, "
        f"mean_bin={outputs['mean_bin'].height}, "
        f"fitted_mean_bin={outputs['fitted_mean_bin'].height}"
    )
    if synonymous is not None:
        reads_per_copy = _reads_per_copy(reads) if synonymous.weighting == "reads" else None
        outputs.update(_combine_synonymous(outputs["per_clonotype"], synonymous, conc_bounds, reads_per_copy))
    return outputs


def _reads_per_copy(reads: pl.DataFrame) -> dict[str, float]:
    """Each nt variant's total reads over the samples of the titration, after the antigen filter."""
    totals = reads.group_by(pl.col(COL_CLONOTYPE).cast(pl.Utf8)).agg(pl.col(COL_READS).fill_null(0).sum())
    return dict(zip(totals[COL_CLONOTYPE].to_list(), totals[COL_READS].cast(pl.Float64).to_list()))


def _combine_synonymous(
    per_clonotype: pl.DataFrame,
    synonymous: SynonymousInputs,
    conc_bounds: tuple[float | None, float | None],
    reads_per_copy: dict[str, float] | None,
) -> dict:
    """The amino-acid outputs: each nt variant's own fit, combined per amino-acid variant."""
    lo, hi = conc_bounds
    per_variant, per_nt_variant, summary = combine_synonymous(
        per_clonotype,
        synonymous,
        log10_lo=math.log10(lo) if lo is not None else None,
        log10_hi=math.log10(hi) if hi is not None else None,
        reads=reads_per_copy,
    )
    weighting = "each nt variant weighted by its total reads" if reads_per_copy is not None else "unweighted"
    log(f"Synonymous combination: mean of log10 Kd,app, {weighting}")
    log(
        f"Synonymous variants: {summary['ntVariants']} nucleotide variants -> {summary['aaVariants']} "
        f"amino-acid variants; left out: {summary['unmappedNtVariants']} with no amino-acid variant, "
        f"{summary['ambiguousNtVariants']} with several"
    )
    if summary["parentMode"]:
        log(
            f"Parent delta: {summary['parents']} parents, {len(summary['unidentifiableParents'])} not identifiable; "
            f"{summary['multiParentVariants']} variants with several parents get no delta"
        )
    else:
        log("Parent delta: off (no amino-acid mutation count and parent link)")
    if "stopAaVariants" in summary:
        log(f"Stop codons: {summary['stopAaVariants']} amino-acid variants not productive; they get no Kd,app")
    else:
        log("Stop codons: off (no amino-acid productive column)")
    outputs = {"per_variant": per_variant, "synonymous_summary": summary}
    if per_nt_variant is not None:
        outputs["per_nt_variant"] = per_nt_variant
    return outputs
