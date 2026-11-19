"""The downloadable copies of the output tables: each variant named by its Variant Id and sequence,
each column headed by the label the results table shows.

The tables the workflow imports keep the hash keys and the column ids, which the PColumn import
reads. The copies offered for download replace every key column with the variant's Variant Id and
sequence, and rename every column the workflow gives a label.
"""

from __future__ import annotations

import polars as pl

from constants import COL_CLONOTYPE
from synonymous import COL_AA, COL_NT

NAME_KEY = "key"
NAME_ID = "variantId"
NAME_SEQUENCE = "sequence"


def read_names(path: str | None) -> pl.DataFrame | None:
    """A workflow-built TSV [key, variantId?, sequence?]. Undefined when the run sent none."""
    if path is None:
        return None
    frame = pl.read_csv(path, separator="\t", infer_schema=False)
    return frame.select([c for c in (NAME_KEY, NAME_ID, NAME_SEQUENCE) if c in frame.columns])


def _named(frame: pl.DataFrame, key: str, names: pl.DataFrame | None, prefix: str) -> pl.DataFrame:
    """`key` replaced in place by the Variant Id and sequence. A key with no Variant Id keeps itself."""
    id_col = f"{prefix}VariantId" if prefix else NAME_ID
    seq_col = f"{prefix}Sequence" if prefix else NAME_SEQUENCE
    keys = pl.col(key).cast(pl.Utf8)
    ids: pl.Expr = keys
    sequences: pl.Expr = pl.lit(None, dtype=pl.Utf8)
    if names is not None:
        if NAME_ID in names.columns:
            mapping = dict(zip(names[NAME_KEY].to_list(), names[NAME_ID].to_list(), strict=True))
            ids = pl.coalesce(keys.replace_strict(mapping, default=None, return_dtype=pl.Utf8), keys)
        if NAME_SEQUENCE in names.columns:
            mapping = dict(zip(names[NAME_KEY].to_list(), names[NAME_SEQUENCE].to_list(), strict=True))
            sequences = keys.replace_strict(mapping, default=None, return_dtype=pl.Utf8)
    at = frame.columns.index(key)
    named = frame.with_columns(ids.alias(id_col), sequences.alias(seq_col))
    order = [*frame.columns[:at], id_col, seq_col, *frame.columns[at + 1 :]]
    return named.select(order)


def name_variants(
    frame: pl.DataFrame, variant_names: pl.DataFrame | None, aa_names: pl.DataFrame | None
) -> pl.DataFrame:
    """Every key column of `frame` replaced by its variant's Variant Id and sequence.

    `variant_names` names the abundance's variants (`clonotypeKey`, and `ntKey` in synonymous mode);
    `aa_names` the amino-acid variants (`aaKey`). A table keyed by both grains prefixes the columns
    "nt" and "aa".
    """
    both = COL_AA in frame.columns and (COL_CLONOTYPE in frame.columns or COL_NT in frame.columns)
    for key, names, prefix in (
        (COL_CLONOTYPE, variant_names, "nt" if both else ""),
        (COL_NT, variant_names, "nt" if both else ""),
        (COL_AA, aa_names, "aa" if both else ""),
    ):
        if key in frame.columns:
            frame = _named(frame, key, names, prefix)
    return frame


def labelled(frame: pl.DataFrame, labels: dict[str, str] | None) -> pl.DataFrame:
    """Each column renamed to its label; a column with no label keeps its name.

    Two columns that would share a header keep them apart by the later one's own name in brackets.
    """
    if not labels:
        return frame
    used: set[str] = set()
    renames: dict[str, str] = {}
    for name in frame.columns:
        header = labels.get(name, name)
        if header in used:
            header = f"{header} ({name})"
        used.add(header)
        renames[name] = header
    return frame.rename(renames)
