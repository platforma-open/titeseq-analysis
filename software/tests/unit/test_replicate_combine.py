"""The replicate combination: what enters the mean, what the error is, and how the mean is rated.

LO and HI are log10 of the lowest and highest concentrations of a titration from 1e-10 to 1e-6.
"""

from __future__ import annotations

import math
import statistics

import polars as pl
import pytest

from replicate_combine import (
    ReplicateValue,
    build_per_replicate_frame,
    collect_values,
    combine_replicates,
    merge_all,
    replicate_value,
    require_synonymous_majority,
)

LO = -10.0
HI = -6.0


def _row(kd: float | None, *, cls: str = "Good", out_of_range: bool = False, se: float | None = None, r2=None):
    return {"kd": kd, "affinityClass": cls, "kdOutOfRange": out_of_range, "kdLog10Se": se, "r2": r2}


def _value(label: str, log10_kd: float, *, se: float | None = None, state="measured", cls="Good"):
    return ReplicateValue(state, log10_kd, se, cls, replicate=label)


def _failed(label: str):
    return ReplicateValue("failed", None, None, "Failed", replicate=label)


class TestValueAndPinning:
    def test_a_measured_replicate_enters_at_its_own_value(self):
        v = replicate_value("A", _row(1e-8, se=0.2), log10_lo=LO, log10_hi=HI)
        assert v.state == "measured"
        assert v.log10_kd == pytest.approx(-8.0)
        assert v.log10_kd_se == 0.2

    def test_a_replicate_past_the_low_wall_enters_pinned_at_that_wall(self):
        v = replicate_value("A", _row(1e-14, out_of_range=True), log10_lo=LO, log10_hi=HI)
        assert (v.state, v.log10_kd) == ("below", LO)

    def test_a_replicate_past_the_high_wall_enters_pinned_at_that_wall(self):
        v = replicate_value("A", _row(1e-2, out_of_range=True), log10_lo=LO, log10_hi=HI)
        assert (v.state, v.log10_kd) == ("above", HI)

    def test_a_pinned_replicate_pulls_the_mean_towards_its_wall(self):
        unpinned = combine_replicates([_value("A", -8.0), _value("B", -8.0)])
        pinned = combine_replicates([_value("A", -8.0), _value("B", LO, state="below")])
        assert pinned.log10_kd < unpinned.log10_kd
        assert pinned.log10_kd == pytest.approx((-8.0 + LO) / 2)

    def test_a_failed_replicate_is_left_out_and_counted(self):
        m = combine_replicates([_value("A", -8.0), _failed("B")])
        assert m.log10_kd == pytest.approx(-8.0)
        assert m.counts["failed"] == 1
        assert m.counts["measured"] == 1

    def test_every_replicate_failed_is_named(self):
        m = combine_replicates([_failed("A"), _failed("B")])
        assert m.log10_kd is None
        assert m.classification.failure_reason == "all_replicates_failed"

    def test_a_non_finite_error_reads_as_no_error(self):
        v = replicate_value("A", _row(1e-8, se=float("nan")), log10_lo=LO, log10_hi=HI)
        assert v.log10_kd_se is None
        assert v.state == "measured"

    def test_a_replicate_value_carries_its_label_r2_and_failure_reason(self):
        row = {**_row(1e-8, se=0.2, r2=0.93), "fitFailureReason": None}
        v = replicate_value("A", row, log10_lo=LO, log10_hi=HI)
        assert (v.replicate, v.r2, v.failure_reason) == ("A", 0.93, None)


class TestError:
    POINTS = [-8.0, -7.5, -8.6]

    def test_the_error_is_the_sample_sd_over_root_n(self):
        m = combine_replicates([_value(str(i), p) for i, p in enumerate(self.POINTS)])
        assert m.log10_kd_se == pytest.approx(statistics.stdev(self.POINTS) / math.sqrt(3))

    def test_one_replicate_keeps_its_own_fit_error(self):
        m = combine_replicates([_value("A", -8.0, se=0.31)])
        assert m.log10_kd_se == 0.31
        assert m.log10_kd_sd is None

    def test_one_replicate_with_no_fit_error_reports_none_and_keeps_its_own_class(self):
        m = combine_replicates([_value("A", -8.0, se=None, cls="Partial")])
        assert m.log10_kd_se is None
        assert m.classification.affinity_class == "Partial"

    def test_one_replicate_with_a_wide_fit_error_rates_partial(self):
        m = combine_replicates([_value("A", -8.0, se=0.61)])
        assert m.log10_kd_se == pytest.approx(0.61)
        assert m.classification.affinity_class == "Partial"

    def test_the_sd_column_is_the_sample_sd(self):
        values = [_value(str(i), p) for i, p in enumerate(self.POINTS)]
        assert combine_replicates(values).log10_kd_sd == pytest.approx(statistics.stdev(self.POINTS))

    def test_a_wide_spread_rates_partial(self):
        m = combine_replicates([_value("A", -9.0), _value("B", -6.0)])
        assert m.log10_kd_se > 0.5
        assert m.classification.affinity_class == "Partial"


class TestFitErrorFloor:
    """A combined SE is never smaller than what the replicates' own fit errors support."""

    def test_close_replicates_with_wide_fits_are_floored_by_their_fit_error(self):
        m = combine_replicates([_value("A", -8.00, se=0.9), _value("B", -8.05, se=0.9)])
        assert m.log10_kd_se == pytest.approx(0.9 / math.sqrt(2))
        assert m.classification.affinity_class == "Partial"

    def test_the_spread_stands_when_it_exceeds_the_fit_error_term(self):
        points = [-9.0, -7.0]
        m = combine_replicates([_value("A", p, se=0.1) for p in points])
        assert m.log10_kd_se == pytest.approx(statistics.stdev(points) / math.sqrt(2))

    def test_the_fit_error_term_is_the_rms_over_root_n_counting_every_entered_replicate(self):
        values = [_value("A", -6.1, se=0.3), _value("B", -6.1, se=0.4), _value("C", HI, state="above", cls="Partial")]
        m = combine_replicates(values)
        assert m.log10_kd_se == pytest.approx(math.sqrt((0.3**2 + 0.4**2) / 2) / math.sqrt(3))

    def test_the_same_members_rate_alike_as_replicates_and_as_synonymous_variants(self):
        from classify import Member
        from synonymous import combine

        replicates = combine_replicates([_value("A", -8.00, se=0.9), _value("B", -8.05, se=0.9)])
        synonymous = combine([Member("measured", -8.00, 0.9, "Good"), Member("measured", -8.05, 0.9, "Good")])
        assert replicates.log10_kd_se == pytest.approx(synonymous.log10_kd_se)
        assert replicates.classification == synonymous.classification


class TestNoFitError:
    """With no fit SE (the mean bin fit), the spread is the only evidence: Good also needs more than
    half the entered replicates rated Good."""

    def test_agreeing_replicates_rated_partial_stay_partial(self):
        m = combine_replicates([_value("A", -8.0, cls="Partial"), _value("B", -8.02, cls="Partial")])
        assert m.log10_kd_se <= 0.5
        assert m.classification.affinity_class == "Partial"

    def test_half_rated_good_is_not_a_majority(self):
        m = combine_replicates([_value("A", -8.0, cls="Good"), _value("B", -8.02, cls="Partial")])
        assert m.classification.affinity_class == "Partial"

    def test_most_rated_good_rates_good(self):
        values = [_value("A", -8.0), _value("B", -8.02), _value("C", -8.01, cls="Partial")]
        assert combine_replicates(values).classification.affinity_class == "Good"

    def test_one_replicate_keeps_its_own_class(self):
        assert combine_replicates([_value("A", -8.0, cls="Good")]).classification.affinity_class == "Good"


class TestCensoring:
    def test_every_replicate_censored_low_reports_that_wall_with_no_error(self):
        m = combine_replicates([_value("A", LO, state="below"), _value("B", LO, state="below")])
        assert (m.log10_kd, m.log10_kd_se, m.out_of_range) == (LO, None, True)
        assert m.classification.affinity_class == "Below range"

    def test_every_replicate_censored_high_reports_that_wall(self):
        m = combine_replicates([_value("A", HI, state="above")])
        assert (m.log10_kd, m.out_of_range) == (HI, True)

    def test_opposite_walls_with_no_measured_replicate_fail(self):
        m = combine_replicates([_value("A", LO, state="below"), _value("B", HI, state="above")])
        assert m.log10_kd is None
        assert m.classification.affinity_class == "Failed"
        assert m.classification.failure_reason == "conflicting_censoring"

    def test_a_majority_on_one_wall_sets_the_bound(self):
        m = combine_replicates(
            [_value("A", LO, state="below"), _value("B", LO, state="below"), _value("C", HI, state="above")]
        )
        assert (m.log10_kd, m.out_of_range) == (LO, True)

    def test_an_even_split_between_the_walls_fails(self):
        values = [_value(r, LO, state="below") for r in ("A", "B")] + [_value(r, HI, state="above") for r in "CD"]
        assert combine_replicates(values).classification.failure_reason == "conflicting_censoring"

    def test_one_measured_replicate_beside_opposite_walls_still_takes_the_mean(self):
        m = combine_replicates([_value("A", -8.0), _value("B", LO, state="below"), _value("C", HI, state="above")])
        assert m.log10_kd == pytest.approx((-8.0 + LO + HI) / 3)
        assert m.out_of_range is False


class TestFrames:
    def _collected(self):
        import polars as pl

        frames = [
            (
                label,
                pl.DataFrame(
                    {
                        "clonotypeKey": ["v1", "v2"],
                        "kd": [kd, None],
                        "affinityClass": ["Good", "Failed"],
                        "kdOutOfRange": [False, None],
                        "kdLog10Se": [0.2, None],
                        "r2": [0.95, None],
                        "hillCoefficient": [1.0, None],
                    }
                ),
            )
            for label, kd in (("A", 1e-8), ("B", 2e-8))
        ]
        return collect_values(frames, "clonotypeKey", log10_lo=LO, log10_hi=HI)

    def test_the_per_replicate_frame_holds_one_row_a_replicate(self):
        frame = build_per_replicate_frame(self._collected(), "clonotypeKey")
        assert frame.columns == [
            "clonotypeKey", "replicate", "replicateReads", "kd", "kdNegLog10", "kdLog10Se", "r2",
            "affinityClass", "fitFailureReason",
        ]  # fmt: skip
        assert frame.height == 4
        assert frame.filter(affinityClass="Failed").height == 2

    def test_the_per_replicate_frame_reports_kd_in_the_published_form(self):
        frame = build_per_replicate_frame(self._collected(), "clonotypeKey")
        measured = frame.filter(pl.col("affinityClass") != "Failed").sort("replicate")
        assert measured["kd"].to_list() == pytest.approx([1e-8, 2e-8])
        assert measured["kdNegLog10"].to_list() == pytest.approx([8.0, 7.69897], abs=1e-5)

    def test_the_counts_add_up_to_the_replicates_seen(self):
        merged = merge_all(self._collected())
        assert merged["v1"].counts == {"measured": 2, "below": 0, "above": 0, "failed": 0}
        assert merged["v2"].counts["failed"] == 2


class TestBoundClassAndRangeSide:
    def test_every_replicate_at_the_high_wall_combines_as_above_range(self):
        values = [_value(r, HI, state="above", cls="Above range") for r in ("A", "B")]
        merged = combine_replicates(values)
        assert merged.classification.affinity_class == "Above range"
        assert merged.range_side == "Above range"

    def test_the_wall_decides_the_class_whatever_each_replicate_was_rated(self):
        values = [
            _value("A", HI, state="above", cls="Above range"),
            _value("B", HI, state="above", cls="Partial"),
        ]
        assert combine_replicates(values).classification.affinity_class == "Above range"

    def test_a_low_wall_bound_is_below_range(self):
        merged = combine_replicates([_value(r, LO, state="below", cls="Below range") for r in ("A", "B")])
        assert (merged.classification.affinity_class, merged.range_side) == ("Below range", "Below range")

    def test_a_mean_of_measured_values_is_in_range(self):
        merged = combine_replicates([_value("A", -8.0, se=0.1), _value("B", -8.2, se=0.1)])
        assert merged.range_side == "In range"

    def test_no_value_has_no_side(self):
        merged = combine_replicates([_value("A", -8.0, state="failed", cls="Failed")])
        assert merged.range_side is None

    def test_the_mean_bin_fit_emits_no_wall_class_for_a_bound(self):
        # "hill_gate_index" declares Good, Partial and Failed only: a wall class would fall outside
        # the column's declared values and a downstream filter would drop the variant unseen.
        values = [_value(r, HI, state="above", cls="Good") for r in ("A", "B")]
        merged = combine_replicates(values, wall_classes=False)
        assert merged.classification.affinity_class == "Good"
        assert merged.range_side == "Above range"

    def test_under_the_mean_bin_fit_replicates_that_disagree_give_partial(self):
        values = [_value("A", HI, state="above", cls="Good"), _value("B", HI, state="above", cls="Partial")]
        assert combine_replicates(values, wall_classes=False).classification.affinity_class == "Partial"


class TestPartlyCensored:
    def test_replicates_at_a_wall_cap_the_class_at_partial(self):
        # Two flat replicates enter at the top wall and agree exactly; the spread is small, but
        # it measures the clipping, not the error.
        values = [
            _value("rep1", -5.0, state="above", cls="Above range"),
            _value("rep2", -5.91, se=0.2),
            _value("rep3", -5.0, state="above", cls="Above range"),
        ]
        m = combine_replicates(values)
        assert m.log10_kd == pytest.approx((-5.0 - 5.91 - 5.0) / 3)
        assert m.log10_kd_se <= 0.5
        assert m.classification.affinity_class == "Partial"

    def test_measured_replicates_alone_still_rate_good(self):
        m = combine_replicates([_value("A", -8.0), _value("B", -8.1), _value("C", -7.9)])
        assert m.classification.affinity_class == "Good"

    def test_a_measured_majority_rates_good_with_one_replicate_at_a_wall(self):
        values = [
            _value("rep1", -9.8),
            _value("rep2", -9.9),
            _value("rep3", -10.0, state="below", cls="Below range"),
        ]
        m = combine_replicates(values)
        assert m.log10_kd == pytest.approx((-9.8 - 9.9 - 10.0) / 3)
        assert m.classification.affinity_class == "Good"

    def test_half_measured_is_not_a_majority(self):
        values = [
            _value("rep1", -9.8),
            _value("rep2", -9.9),
            _value("rep3", -10.0, state="below", cls="Below range"),
            _value("rep4", -10.0, state="below", cls="Below range"),
        ]
        m = combine_replicates(values)
        assert m.log10_kd_se <= 0.5
        assert m.classification.affinity_class == "Partial"


def _aa_row(log10_kd, se, *, measured, below=0, above=0, cls="Partial"):
    """One replicate's amino-acid value, as the synonymous combination writes it."""
    return {
        "kd": 10.0**log10_kd,
        "affinityClass": cls,
        "kdOutOfRange": False,
        "kdLog10Se": se,
        "ntMeasured": measured,
        "ntCensoredBelow": below,
        "ntCensoredAbove": above,
    }


class TestSynonymousMajority:
    """Starr order: the synonymous variants combine inside each replicate, then the replicates combine."""

    def _merged(self, rows):
        frames = [(f"rep{i}", pl.DataFrame([{"aaKey": "A", **row}])) for i, row in enumerate(rows)]
        collected = collect_values(frames, "aaKey", log10_lo=LO, log10_hi=HI)
        return require_synonymous_majority(collected, merge_all(collected))["A"]

    def test_replicates_built_mostly_from_synonymous_variants_at_a_wall_do_not_rate_good(self):
        # Each replicate: one synonymous variant measured near -8, two at the low wall. Each replicate
        # value is measured, and the replicates agree on the walls, so the SE alone would say Good.
        rows = [_aa_row(x, 0.6, measured=1, below=2) for x in (-9.33, -9.36, -9.30)]
        m = self._merged(rows)
        assert m.log10_kd_se <= 0.5
        assert m.classification.affinity_class == "Partial"

    def test_a_synonymous_majority_over_the_replicates_still_rates_good(self):
        rows = [_aa_row(x, 0.2, measured=2, below=1) for x in (-8.65, -8.70, -8.60)]
        assert self._merged(rows).classification.affinity_class == "Good"

    def test_half_the_synonymous_variants_measured_is_not_a_majority(self):
        rows = [_aa_row(-9.0, 0.2, measured=1, below=1), _aa_row(-9.05, 0.2, measured=1, above=1)]
        assert self._merged(rows).classification.affinity_class == "Partial"

    def test_a_bound_keeps_its_class_under_the_mean_bin_fit(self):
        # hill_gate_index rates a fit on R² and n, so synonymous variants past the top wall can be
        # Good. Every replicate's amino-acid value is that bound, as it is without replicates. The
        # majority rule guards a mean's SE, and a bound has none.
        rows = [{**_aa_row(HI, None, measured=0, above=2, cls="Good"), "kdOutOfRange": True} for _ in range(3)]
        frames = [(f"rep{i}", pl.DataFrame([{"aaKey": "A", **row}])) for i, row in enumerate(rows)]
        collected = collect_values(frames, "aaKey", log10_lo=LO, log10_hi=HI)
        m = require_synonymous_majority(collected, merge_all(collected, wall_classes=False))["A"]
        assert m.out_of_range
        assert m.classification.affinity_class == "Good"

    def _merged_spread_only(self, rows):
        frames = [(f"rep{i}", pl.DataFrame([{"aaKey": "A", **row}])) for i, row in enumerate(rows)]
        collected = collect_values(frames, "aaKey", log10_lo=LO, log10_hi=HI)
        return require_synonymous_majority(collected, merge_all(collected, wall_classes=False, fits_report_se=False))[
            "A"
        ]

    def test_under_the_mean_bin_fit_partial_replicates_stay_partial(self):
        # Each replicate's synonymous variants are Partial fits that agree, so its amino-acid value is
        # Partial with a small SE. That SE is a spread, not a fit SE, so the replicates' mean still
        # needs a Good majority.
        rows = [_aa_row(x, 0.01, measured=2, cls="Partial") for x in (-8.0, -8.02)]
        m = self._merged_spread_only(rows)
        assert m.log10_kd_se <= 0.5
        assert m.classification.affinity_class == "Partial"

    def test_under_the_mean_bin_fit_one_partial_replicate_stays_partial(self):
        m = self._merged_spread_only([_aa_row(-8.0, 0.01, measured=2, cls="Partial")])
        assert m.log10_kd_se == 0.01
        assert m.classification.affinity_class == "Partial"

    def test_under_the_mean_bin_fit_good_replicates_rate_good(self):
        rows = [_aa_row(x, 0.01, measured=2, cls="Good") for x in (-8.0, -8.02, -8.01)]
        assert self._merged_spread_only(rows).classification.affinity_class == "Good"

    def test_synonymous_variants_of_a_failed_replicate_are_not_counted(self):
        # The failed replicate's censored synonymous variants never entered the mean.
        failed = {**_aa_row(-8.0, None, measured=0, below=3), "kd": None, "affinityClass": "Failed"}
        rows = [_aa_row(-8.0, 0.2, measured=2), _aa_row(-8.05, 0.2, measured=2), failed]
        assert self._merged(rows).classification.affinity_class == "Good"
