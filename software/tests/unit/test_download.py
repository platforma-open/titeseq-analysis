"""The downloadable copies name each variant by its Variant Id and sequence, not its hash key."""

from __future__ import annotations

import polars as pl
from polars.testing import assert_frame_equal

from download import labelled, name_variants
from synonymous import nt_aa_variants

NAMES = pl.DataFrame({"key": ["h1", "h2"], "variantId": ["V-AAAAA", "V-BBBBB"], "sequence": ["CARW", "CARY"]})
AA_NAMES = pl.DataFrame({"key": ["a1"], "variantId": ["V-CCCCC"], "sequence": ["CAR"]})


def test_the_key_is_replaced_in_place_by_id_and_sequence():
    frame = pl.DataFrame({"clonotypeKey": ["h2", "h1"], "kd": [1e-8, 2e-8]})
    expected = pl.DataFrame({"variantId": ["V-BBBBB", "V-AAAAA"], "sequence": ["CARY", "CARW"], "kd": [1e-8, 2e-8]})
    assert_frame_equal(name_variants(frame, NAMES, None), expected)


def test_amino_acid_tables_use_the_amino_acid_names():
    frame = pl.DataFrame({"aaKey": ["a1"], "replicate": ["r1"]})
    expected = pl.DataFrame({"variantId": ["V-CCCCC"], "sequence": ["CAR"], "replicate": ["r1"]})
    assert_frame_equal(name_variants(frame, NAMES, AA_NAMES), expected)


def test_a_table_keyed_by_both_grains_prefixes_them():
    frame = pl.DataFrame({"ntKey": ["h1"], "aaKey": ["a1"], "synonymous": [True]})
    out = name_variants(frame, NAMES, AA_NAMES)
    assert out.columns == ["ntVariantId", "ntSequence", "aaVariantId", "aaSequence", "synonymous"]
    assert out.row(0) == ("V-AAAAA", "CARW", "V-CCCCC", "CAR", True)


def test_a_variant_with_no_name_keeps_its_key_and_an_empty_sequence():
    frame = pl.DataFrame({"clonotypeKey": ["h3"]})
    out = name_variants(frame, NAMES, None)
    assert out.row(0) == ("h3", None)


def test_no_sequence_column_upstream_leaves_the_sequence_empty():
    frame = pl.DataFrame({"clonotypeKey": ["h1"]})
    out = name_variants(frame, NAMES.drop("sequence"), None)
    assert out.columns == ["variantId", "sequence"]
    assert out.row(0) == ("V-AAAAA", None)


def test_no_names_at_all_keeps_the_key_as_the_id():
    frame = pl.DataFrame({"clonotypeKey": ["h1"], "kd": [1.0]})
    assert name_variants(frame, None, None).row(0) == ("h1", None, 1.0)


def test_a_table_with_no_key_is_unchanged():
    frame = pl.DataFrame({"concentration": [1.0]})
    assert_frame_equal(name_variants(frame, NAMES, AA_NAMES), frame)


def test_the_nt_to_aa_table_names_both_grains():
    frame = nt_aa_variants(pl.DataFrame({"ntKey": ["h1"], "aaKey": ["a1"]}))
    out = name_variants(frame, NAMES, AA_NAMES)
    assert out.columns == ["ntVariantId", "ntSequence", "aaVariantId", "aaSequence"]
    assert out.row(0) == ("V-AAAAA", "CARW", "V-CCCCC", "CAR")


def test_columns_take_their_labels_and_unlabelled_ones_keep_their_names():
    frame = pl.DataFrame({"variantId": ["V-A"], "kd": [1.0], "ntVariants": [2]})
    out = labelled(frame, {"variantId": "Variant Id", "kd": "Kd,app"})
    assert out.columns == ["Variant Id", "Kd,app", "ntVariants"]


def test_two_columns_with_one_label_stay_apart():
    frame = pl.DataFrame({"kdLog10Se": [0.1], "kdLog10Sd": [0.2]})
    assert labelled(frame, {"kdLog10Se": "Kd,app spread", "kdLog10Sd": "Kd,app spread"}).columns == [
        "Kd,app spread",
        "Kd,app spread (kdLog10Sd)",
    ]


def test_no_labels_leave_the_frame_as_it_is():
    frame = pl.DataFrame({"kd": [1.0]})
    assert labelled(frame, None).columns == ["kd"]
