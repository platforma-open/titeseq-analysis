"""CLI entry point — parses args, reads the input table, delegates to pipeline.run, writes outputs."""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

import polars as pl

from constants import DEFAULT_PARAMS, FitParams
from io_layer import read_reads_table
from log import log
from pipeline import run


def _write_frame(frame: pl.DataFrame, path: str) -> None:
    p = Path(path)
    if p.suffix in (".parquet", ".pq"):
        frame.write_parquet(p)
    else:
        frame.write_csv(p, separator="\t")


def fit_summary(per_clonotype: pl.DataFrame) -> dict[str, int]:
    return {
        "clonotypes": per_clonotype.height,
        "failed": per_clonotype.filter(pl.col("affinityClass") == "Failed").height,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="fit-curves")
    parser.add_argument("--reads", required=True, help="path to reads table (parquet or tsv)")
    parser.add_argument("--out-per-clonotype", required=True)
    parser.add_argument("--out-mean-bin", required=True)
    parser.add_argument("--out-fitted-mean-bin", required=True)
    parser.add_argument("--out-concentration-value", required=True)
    parser.add_argument(
        "--out-fit-summary",
        default=None,
        help="Optional JSON with the clonotype count and how many are Failed.",
    )
    parser.add_argument("--params", default=None, help="path to params JSON (optional)")
    parser.add_argument("--target-antigen", default=None)
    parser.add_argument("--antigen-column-ref", default=None)
    parser.add_argument(
        "--sort-fraction-column",
        default=None,
        help="Column name in reads.tsv carrying sort fraction C_bc/C_c per sample. Absent ⇒ no FACS correction.",
    )
    parser.add_argument(
        "--cell-count-column",
        default=None,
        help="Column name in reads.tsv carrying the cells sorted into each sample (C_bc). "
        "Preferred over --sort-fraction-column: also allows rescaling bins with more reads than cells.",
    )
    parser.add_argument(
        "--allow-computed-baseline",
        action="store_true",
        help="Compute the 0 M baseline when none is given. Only correct when every clonotype is a "
        "designed library member: on a de-novo clonotype set the sequencing-error clonotypes each "
        "carry equal weight and inflate it several-fold. Prefer setting baseline_fluorescence.",
    )
    args = parser.parse_args(argv)

    log(f"Reading input table: {args.reads}")
    reads = read_reads_table(args.reads)
    log(f"Loaded {reads.height} rows, {len(reads.columns)} columns")
    params = DEFAULT_PARAMS
    if args.params:
        with open(args.params) as f:
            params = FitParams(**json.load(f))
        log(f"Loaded custom params from {args.params}")
    if args.allow_computed_baseline:
        params = dataclasses.replace(params, allow_computed_baseline=True)

    outputs = run(
        reads,
        params=params,
        target_antigen=args.target_antigen,
        antigen_column_ref=args.antigen_column_ref,
        sort_fraction_col=args.sort_fraction_column,
        cell_count_col=args.cell_count_column,
    )
    log("Writing output frames")
    _write_frame(outputs["per_clonotype"], args.out_per_clonotype)
    _write_frame(outputs["mean_bin"], args.out_mean_bin)
    _write_frame(outputs["fitted_mean_bin"], args.out_fitted_mean_bin)
    _write_frame(outputs["concentration_value"], args.out_concentration_value)
    if args.out_fit_summary:
        Path(args.out_fit_summary).write_text(json.dumps(fit_summary(outputs["per_clonotype"])))
    log("CLI done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
