"""Session fixtures for the Adams 2016 golden tests (data loading lives in golden_data.py)."""

from __future__ import annotations

import io

import polars as pl
import pytest
from golden_data import GoldenReplicate, build_replicate, fetch
from paper_adams2016 import REPLICATES


@pytest.fixture(scope="session")
def adams_replicates() -> list[GoldenReplicate]:
    return [build_replicate(rep) for rep in REPLICATES]


@pytest.fixture(scope="session")
def adams_clone_keys() -> dict[str, str]:
    """Table I clone name -> clonotypeKey (CDR1H_CDR3H nucleotides), from the authors' flow data."""
    # The file uses classic-Mac "\r" line endings.
    text = fetch("titration_curves.csv").read_text().replace("\r\n", "\n").replace("\r", "\n")
    curves = pl.read_csv(io.StringIO(text), infer_schema=False)
    curves = curves.rename({c: c.strip() for c in curves.columns})
    return {
        str(row["name"]).strip(): f"{str(row['CDR1H']).strip()}_{str(row['CDR3H']).strip()}"
        for row in curves.unique(subset="name", keep="first").iter_rows(named=True)
    }
