"""CLI entry point — orchestrates the full fitting pipeline.

Pipeline:
  1. Load reads table.
  2. Validate R1-R5 and apply antigen filter (R4).
  3. Compute signal (R7 bin mode or R7b no-bin mode).
  4. Apply R8 weights + R9 floor filter.
  5. Compute global baseline B from c=0 points (R6).
  6. Detect hook effect per clonotype (R9b).
  7. Fit Hill equation per clonotype (R10) + weighted R² (R11).
  8. Classify (R12), flag kdOutOfRange (R14b).
  9. Write outputs: per-clonotype frame + meanBin + fittedMeanBin (R13, R14).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

from baseline import compute_global_baseline
from classify import classify
from constants import (
    COL_BIN,
    COL_CLONOTYPE,
    COL_CONC_STR,
    COL_CONC_VAL,
    COL_ANTIGEN,
    DEFAULT_PARAMS,
    FitParams,
)
from floor_filter import apply_floor_and_weights, classify_insufficient, split_c0, WEIGHT
from hill_fit import fit_one_clonotype
from hook_effect import detect_hook_effect
from io_layer import (
    apply_antigen_filter,
    canonicalize_concentration,
    max_bin_label,
    read_reads_table,
    validate_antigen_filter,
    validate_bin_column,
    validate_concentration_column,
    validate_reads_schema,
    validate_sample_metadata_uniqueness,
)
from normalization import SIGNAL, normalize
from output_build import (
    FITTED_MEAN_BIN_SCHEMA,
    PER_CLONOTYPE_SCHEMA,
    build_mean_bin_frame,
    flag_kd_out_of_range,
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
    warnings += validate_antigen_filter(reads, antigen_column_ref, target_antigen)
    validate_sample_metadata_uniqueness(reads, has_bin=has_bin, has_antigen=has_antigen)
    for w in warnings:
        print(f"WARN: {w}", file=sys.stderr)


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
    """Iterate clonotypes once; return (per_clonotype_frame, fitted_mean_bin_frame).

    Partitions fit_points once (O(N)) instead of filtering per clonotype inside
    the loop (O(N * C)). Clonotypes failing the pre-fit gates (insufficient_*,
    non_monotonic_signal) skip partition lookup entirely.

    Output frames are built from pre-allocated numpy arrays and fitted-row chunks
    — avoids the C × 7 Python dict allocations (per-clonotype) + Σk_c × 4 dicts
    (fitted rows) that list-of-dicts → DataFrame conversion otherwise incurs.
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

    fitted_keys: list[np.ndarray] = []
    fitted_conc_strs: list[np.ndarray] = []
    fitted_conc_vals: list[np.ndarray] = []
    fitted_y_hats: list[np.ndarray] = []

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
        y = sub[SIGNAL].to_numpy()
        w = sub[WEIGHT].to_numpy()

        fit = fit_one_clonotype(
            x, y, w, baseline_fixed=global_baseline, bin_mode=bin_mode, max_bin_label=max_bin
        )
        cls = classify(fit.r2_w, fit.n, fit.converged, params)

        kd_col[i] = fit.kd
        n_col[i] = fit.n
        r2_col[i] = fit.r2_w
        affinity_col[i] = cls.affinity_class
        reason_col[i] = cls.failure_reason

        if fit.converged and fit.y_hat is not None and cls.affinity_class != "Failed":
            k = x.shape[0]
            fitted_keys.append(np.full(k, clonotype, dtype=object))
            fitted_conc_strs.append(sub[COL_CONC_STR].to_numpy())
            fitted_conc_vals.append(x)
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
                COL_CONC_VAL: np.concatenate(fitted_conc_vals),
                "fittedMeanBin": np.concatenate(fitted_y_hats),
            },
            schema=FITTED_MEAN_BIN_SCHEMA,
        )
    else:
        fitted = pl.DataFrame(
            {k: [] for k in FITTED_MEAN_BIN_SCHEMA}, schema=FITTED_MEAN_BIN_SCHEMA
        )

    return per_clonotype, fitted


def _build_outputs(
    per_clonotype: pl.DataFrame,
    fitted: pl.DataFrame,
    signal_frame: pl.DataFrame,
    reads: pl.DataFrame,
) -> dict[str, pl.DataFrame]:
    """Apply R14b kdOutOfRange flag and assemble the three-frame output dict."""
    min_max = reads.filter(pl.col(COL_CONC_VAL) > 0).select(
        pl.col(COL_CONC_VAL).min().alias("min"),
        pl.col(COL_CONC_VAL).max().alias("max"),
    ).row(0, named=True)
    if min_max["min"] is not None:
        per_clonotype = flag_kd_out_of_range(
            per_clonotype, float(min_max["min"]), float(min_max["max"])
        )

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
) -> dict[str, pl.DataFrame]:
    """Execute the full pipeline on an in-memory reads frame."""
    reads = canonicalize_concentration(reads)
    has_bin = COL_BIN in reads.columns
    has_antigen = COL_ANTIGEN in reads.columns

    _validate_inputs(
        reads, has_bin=has_bin, has_antigen=has_antigen,
        target_antigen=target_antigen, antigen_column_ref=antigen_column_ref,
    )

    reads = apply_antigen_filter(reads, target_antigen)
    mbl = max_bin_label(reads) if has_bin else None

    signal_frame = normalize(reads, bin_mode=has_bin)
    floor_frame = apply_floor_and_weights(signal_frame, params)

    all_clonotypes = signal_frame[COL_CLONOTYPE].unique().to_list()
    insufficient = classify_insufficient(floor_frame, all_clonotypes, params)
    insufficient_map = dict(
        zip(insufficient[COL_CLONOTYPE].to_list(), insufficient["insufficient_reason"].to_list())
    )

    c0_points, fit_points = split_c0(floor_frame)
    global_b = compute_global_baseline(c0_points)

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

    return _build_outputs(per_clonotype, fitted, signal_frame, reads)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="fit-curves")
    parser.add_argument("--reads", required=True, help="path to reads table (parquet or tsv)")
    parser.add_argument("--out-per-clonotype", required=True)
    parser.add_argument("--out-mean-bin", required=True)
    parser.add_argument("--out-fitted-mean-bin", required=True)
    parser.add_argument("--params", default=None, help="path to params JSON (optional)")
    parser.add_argument("--target-antigen", default=None)
    parser.add_argument("--antigen-column-ref", default=None)
    args = parser.parse_args(argv)

    reads = read_reads_table(args.reads)
    params = DEFAULT_PARAMS
    if args.params:
        with open(args.params) as f:
            params = FitParams(**json.load(f))

    outputs = run(
        reads,
        params=params,
        target_antigen=args.target_antigen,
        antigen_column_ref=args.antigen_column_ref,
    )
    _write_frame(outputs["per_clonotype"], args.out_per_clonotype)
    _write_frame(outputs["mean_bin"], args.out_mean_bin)
    _write_frame(outputs["fitted_mean_bin"], args.out_fitted_mean_bin)
    return 0


def _write_frame(frame: pl.DataFrame, path: str) -> None:
    p = Path(path)
    if p.suffix in (".parquet", ".pq"):
        frame.write_parquet(p)
    else:
        frame.write_csv(p, separator="\t")


if __name__ == "__main__":
    sys.exit(main())
