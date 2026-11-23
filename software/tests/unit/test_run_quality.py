import polars as pl

from run_quality import failure_reason_counts, fit_grain, measured, outcome_counts, raw_counts, read_known, summary


def _frame(with_replicate: bool) -> pl.DataFrame:
    rows = {
        "clonotypeKey": ["a", "b", "c", "d", "e", "f"],
        "affinityClass": ["Good", "Good", "Failed", "Failed", "Failed", "Partial"],
        "fitFailureReason": [None, None, "insufficient_points", "insufficient_points", "insufficient_reads", None],
    }
    if with_replicate:
        rows["replicate"] = ["1", "2", "1", "1", "2", "2"]
    return pl.DataFrame(rows)


def test_outcomes_without_replicates_carry_no_replicate():
    out = outcome_counts(_frame(False))
    assert out.columns == ["affinityClass", "count"]
    assert dict(out.iter_rows()) == {"Failed": 3, "Good": 2, "Partial": 1}


def test_outcomes_per_replicate_follow_the_counts_over_all_replicates():
    out = outcome_counts(_frame(True))
    assert out.rows()[:3] == [
        ("All replicates", "Failed", 3),
        ("All replicates", "Good", 2),
        ("All replicates", "Partial", 1),
    ]
    assert out.rows()[3:] == [
        ("1", "Failed", 2),
        ("1", "Good", 1),
        ("2", "Failed", 1),
        ("2", "Good", 1),
        ("2", "Partial", 1),
    ]


def test_failure_reasons_are_labelled_and_cover_failed_fits_only():
    out = failure_reason_counts(_frame(False))
    assert out.rows() == [("Insufficient points", 2), ("Insufficient reads", 1)]


def test_no_failures_gives_an_empty_table():
    frame = _frame(True).filter(pl.col("affinityClass") != "Failed")
    out = failure_reason_counts(frame)
    assert out.columns == ["replicate", "failureReason", "count"]
    assert out.height == 0


def test_fit_grain_prefers_the_lowest_grain():
    a, b, c = pl.DataFrame({"x": [1]}), pl.DataFrame({"x": [2]}), pl.DataFrame({"x": [3]})
    assert fit_grain({"per_clonotype": a}) is a
    assert fit_grain({"per_clonotype": a, "per_replicate": b}) is b
    assert fit_grain({"per_clonotype": a, "per_replicate": b, "per_replicate_nt": c}) is c


def test_without_a_known_list_the_fits_never_attempted_are_left_out():
    kept = measured(_frame(False), None)
    assert dict(outcome_counts(kept).iter_rows()) == {"Good": 2, "Partial": 1}
    assert failure_reason_counts(kept).height == 0


def test_with_a_known_list_only_its_variants_count_whatever_their_outcome():
    kept = measured(_frame(False), {"a", "c", "e", "zz"})
    assert dict(outcome_counts(kept).iter_rows()) == {"Failed": 2, "Good": 1}
    assert failure_reason_counts(kept).rows() == [("Insufficient points", 1), ("Insufficient reads", 1)]


def test_known_variants_are_the_exact_matches(tmp_path):
    path = tmp_path / "known.tsv"
    path.write_text("key\tdistance\na\t0\nb\t1\nc\t0\nd\t\n")
    assert read_known(str(path)) == {"a", "c"}
    assert read_known(None) is None


def test_raw_counts_keep_every_fit():
    assert raw_counts(_frame(False)).rows() == [
        ("Failed", "Insufficient points", 2),
        ("Failed", "Insufficient reads", 1),
        ("Good", None, 2),
        ("Partial", None, 1),
    ]


def test_summary_rows():
    frame = _frame(False)
    s = summary(outcome_counts(frame), failure_reason_counts(frame), raw_counts(frame), None)
    assert {"affinityClass": "Good", "count": 2} in s["outcomes"]
    assert {"failureReason": "Insufficient reads", "count": 1} in s["failureReasons"]
    assert {"affinityClass": "Good", "failureReason": None, "count": 2} in s["raw"]
    assert s["scope"] == "attempted"
    assert s["leftOutReasons"] == ["Insufficient reads", "Insufficient points"]


def test_summary_with_a_known_list_leaves_no_reason_out():
    frame = _frame(False)
    s = summary(outcome_counts(frame), failure_reason_counts(frame), raw_counts(frame), {"a"})
    assert s["scope"] == "known"
    assert s["leftOutReasons"] == []
