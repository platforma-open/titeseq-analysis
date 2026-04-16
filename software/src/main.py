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
from normalization import CLONOTYPE_READS_AT_CONC, SIGNAL, normalize
from output_build import (
    build_fitted_mean_bin_frame,
    build_mean_bin_frame,
    build_per_clonotype_frame,
    flag_kd_out_of_range,
)


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

    validate_reads_schema(reads, has_bin=has_bin, has_antigen=has_antigen)
    validate_concentration_column(reads, has_bin=has_bin)
    if has_bin:
        validate_bin_column(reads)
    validate_antigen_filter(reads, antigen_column_ref, target_antigen)
    validate_sample_metadata_uniqueness(reads, has_bin=has_bin, has_antigen=has_antigen)

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

    per_clonotype_rows: list[dict] = []
    fitted_rows: list[dict] = []

    for clonotype in sorted(all_clonotypes):
        insufficient_reason = insufficient_map.get(clonotype)
        if insufficient_reason is not None:
            per_clonotype_rows.append(
                _failed_row(clonotype, insufficient_reason)
            )
            continue
        if hook_map.get(clonotype, False):
            per_clonotype_rows.append(_failed_row(clonotype, "non_monotonic_signal"))
            continue

        sub = fit_points.filter(pl.col(COL_CLONOTYPE) == clonotype).sort(COL_CONC_VAL)
        x = sub[COL_CONC_VAL].to_numpy()
        y = sub[SIGNAL].to_numpy()
        w = sub[WEIGHT].to_numpy()

        fit = fit_one_clonotype(
            x, y, w, baseline_fixed=global_b, bin_mode=has_bin, max_bin_label=mbl
        )
        cls = classify(fit.r2_w, fit.n, fit.converged, params)
        per_clonotype_rows.append(
            {
                COL_CLONOTYPE: clonotype,
                "kd": fit.kd,
                "hillCoefficient": fit.n,
                "r2": fit.r2_w,
                "affinityClass": cls.affinity_class,
                "fitFailureReason": cls.failure_reason,
                "kdOutOfRange": None,
            }
        )

        if (
            fit.converged
            and fit.y_hat is not None
            and cls.affinity_class != "Failed"
        ):
            conc_strs = sub[COL_CONC_STR].to_list()
            for i, cs in enumerate(conc_strs):
                fitted_rows.append(
                    {
                        COL_CLONOTYPE: clonotype,
                        COL_CONC_STR: cs,
                        COL_CONC_VAL: float(x[i]),
                        "fittedMeanBin": float(fit.y_hat[i]),
                    }
                )

    per_clonotype = build_per_clonotype_frame(per_clonotype_rows)
    # Flag kdOutOfRange from experimental concentration range (exclude c=0).
    non_zero = reads.filter(pl.col(COL_CONC_VAL) > 0)
    if non_zero.height > 0:
        min_c = float(non_zero[COL_CONC_VAL].min())
        max_c = float(non_zero[COL_CONC_VAL].max())
        per_clonotype = flag_kd_out_of_range(per_clonotype, min_c, max_c)

    mean_bin_out = build_mean_bin_frame(signal_frame)
    fitted_mean_bin_out = build_fitted_mean_bin_frame(fitted_rows)

    return {
        "per_clonotype": per_clonotype,
        "mean_bin": mean_bin_out,
        "fitted_mean_bin": fitted_mean_bin_out,
    }


def _failed_row(clonotype: str, reason: str) -> dict:
    return {
        COL_CLONOTYPE: clonotype,
        "kd": None,
        "hillCoefficient": None,
        "r2": None,
        "affinityClass": "Failed",
        "fitFailureReason": reason,
        "kdOutOfRange": None,
    }


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
