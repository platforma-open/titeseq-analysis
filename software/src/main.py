"""CLI entry point — parses args, reads the input table, delegates to pipeline.run, writes outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import polars as pl

from constants import DEFAULT_PARAMS, FitParams
from io_layer import read_reads_table
from pipeline import run


def _write_frame(frame: pl.DataFrame, path: str) -> None:
    p = Path(path)
    if p.suffix in (".parquet", ".pq"):
        frame.write_parquet(p)
    else:
        frame.write_csv(p, separator="\t")


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


if __name__ == "__main__":
    sys.exit(main())
