"""CLI entry point — parses args, reads the input table, delegates to pipeline.run, writes outputs."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import polars as pl

from constants import DEFAULT_PARAMS, FitParams
from io_layer import read_reads_table
from pipeline import run


def _available_cpus() -> int:
    """CPU budget respecting Linux cgroup/cpuset affinity when available.

    os.cpu_count() returns the host count, which over-reports inside containers
    with CPU limits (Docker, K8s). os.sched_getaffinity is the affinity-aware
    alternative — only exposed on Linux/FreeBSD builds.
    """
    getaffinity = getattr(os, "sched_getaffinity", None)
    if getaffinity is not None:
        return len(getaffinity(0))
    return os.cpu_count() or 1


def _resolve_workers(value: str | None) -> int:
    """Map CLI --workers (str or None) to an int. 'auto' -> available CPU budget."""
    if value is None:
        return 1
    if value == "auto":
        return _available_cpus()
    n = int(value)
    if n < 1:
        raise ValueError(f"--workers must be >= 1 or 'auto', got {value!r}")
    return n


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
    parser.add_argument(
        "--workers",
        default=None,
        help="fit worker count; int or 'auto' (os.cpu_count). Defaults to 1 (serial).",
    )
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
        workers=_resolve_workers(args.workers),
    )
    _write_frame(outputs["per_clonotype"], args.out_per_clonotype)
    _write_frame(outputs["mean_bin"], args.out_mean_bin)
    _write_frame(outputs["fitted_mean_bin"], args.out_fitted_mean_bin)
    return 0


if __name__ == "__main__":
    sys.exit(main())
