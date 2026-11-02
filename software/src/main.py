"""CLI entry point — parses args, reads the input table, delegates to pipeline.run, writes outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

from constants import DEFAULT_CENSORED_MAJORITY, DEFAULT_PARAMS, REPLICATE_ERROR_MODELS, FitParams
from io_layer import read_reads_table
from log import log
from pipeline import run
from replicate_combine import Filters, ReplicateInputs


def _write_frame(frame: pl.DataFrame, path: str) -> None:
    p = Path(path)
    if p.suffix in (".parquet", ".pq"):
        frame.write_parquet(p)
    else:
        frame.write_csv(p, separator="\t")


_REPLICATE_ONLY = ("out_per_replicate", "replicate_fit_error_cutoff", "replicate_r2_floor")


def _replicate_inputs(args: argparse.Namespace, parser: argparse.ArgumentParser) -> ReplicateInputs | None:
    if args.replicate_column is None:
        if any(getattr(args, name) is not None for name in _REPLICATE_ONLY):
            parser.error("the replicate outputs and filters need --replicate-column")
        if args.replicate_error_model != REPLICATE_ERROR_MODELS[0]:
            parser.error("--replicate-error-model needs --replicate-column")
        if args.replicate_censored_majority != DEFAULT_CENSORED_MAJORITY:
            parser.error("--replicate-censored-majority needs --replicate-column")
        return None
    if not args.out_per_replicate:
        parser.error("--replicate-column needs --out-per-replicate")
    for name, value in (
        ("--replicate-fit-error-cutoff", args.replicate_fit_error_cutoff),
        ("--replicate-r2-floor", args.replicate_r2_floor),
    ):
        if value is not None and value < 0:
            parser.error(f"{name} must not be negative")
    if not 0.5 <= args.replicate_censored_majority <= 1.0:
        parser.error("--replicate-censored-majority must be in [0.5, 1]")
    return ReplicateInputs(
        column=args.replicate_column,
        error_model=args.replicate_error_model,
        filters=Filters(fit_error_cutoff=args.replicate_fit_error_cutoff, r2_floor=args.replicate_r2_floor),
        censored_majority=args.replicate_censored_majority,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="fit-curves")
    parser.add_argument("--reads", required=True, help="path to reads table (parquet or tsv)")
    parser.add_argument("--out-per-clonotype", required=True)
    parser.add_argument("--out-mean-bin", required=True)
    parser.add_argument("--out-fitted-mean-bin", required=True)
    parser.add_argument("--out-concentration-value", required=True)
    parser.add_argument("--params", default=None, help="path to params JSON (optional)")
    parser.add_argument("--target-antigen", default=None)
    parser.add_argument("--antigen-column-ref", default=None)
    parser.add_argument(
        "--sort-fraction-column",
        default=None,
        help="Column name in reads.tsv carrying sort fraction C_bc/C_c per sample. Absent ⇒ no FACS correction.",
    )
    parser.add_argument(
        "--replicate-column",
        default=None,
        help="Column name in reads.tsv naming each sample's replicate. Absent ⇒ one titration.",
    )
    parser.add_argument(
        "--replicate-error-model",
        choices=REPLICATE_ERROR_MODELS,
        default=REPLICATE_ERROR_MODELS[0],
        help="sd_over_sqrt_n: sample SD / sqrt(n); phillips: population SD / sqrt(n + 1)",
    )
    parser.add_argument(
        "--replicate-fit-error-cutoff",
        type=float,
        default=None,
        help="Decades. Leave out to keep every replicate whatever its fit error.",
    )
    parser.add_argument(
        "--replicate-r2-floor",
        type=float,
        default=None,
        help="Leave out to keep every replicate whatever its R².",
    )
    parser.add_argument(
        "--replicate-censored-majority",
        type=float,
        default=DEFAULT_CENSORED_MAJORITY,
        help="Share of the entered replicates one censoring side must hold to set the bound.",
    )
    parser.add_argument("--out-per-replicate", default=None, help="one row per (variant, replicate)")
    args = parser.parse_args(argv)
    replicates = _replicate_inputs(args, parser)

    log(f"Reading input table: {args.reads}")
    reads = read_reads_table(args.reads)
    log(f"Loaded {reads.height} rows, {len(reads.columns)} columns")
    params = DEFAULT_PARAMS
    if args.params:
        # FitParams rejects an unknown `equation`. An absent one means a pre-existing project on
        # v2.4.0's fit; see FitParams.
        with open(args.params) as f:
            params = FitParams(**json.load(f))
        log(f"Loaded custom params from {args.params}")

    outputs = run(
        reads,
        params=params,
        target_antigen=args.target_antigen,
        antigen_column_ref=args.antigen_column_ref,
        sort_fraction_col=args.sort_fraction_column,
        replicates=replicates,
    )
    log("Writing output frames")
    _write_frame(outputs["per_clonotype"], args.out_per_clonotype)
    _write_frame(outputs["mean_bin"], args.out_mean_bin)
    _write_frame(outputs["fitted_mean_bin"], args.out_fitted_mean_bin)
    _write_frame(outputs["concentration_value"], args.out_concentration_value)
    if args.out_per_replicate is not None:
        _write_frame(outputs["per_replicate"], args.out_per_replicate)
    log("CLI done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
