from __future__ import annotations

import polars as pl

from fit_reads import fit_reads, with_fit_reads


def test_fit_reads_sums_each_key():
    reads = pl.DataFrame({"k": ["a", "a", "b"], "reads": [3, None, 5]})
    assert dict(fit_reads(reads, "k", "reads").iter_rows()) == {"a": 3, "b": 5}


def test_a_key_with_no_reads_gets_zero():
    frame = with_fit_reads(pl.DataFrame({"k": ["a", "b"]}), "k", {"a": 7})
    assert frame["fitReads"].to_list() == [7, 0]
