"""The replicate combination: what enters the mean, what the error is, and what the filters drop.

LO and HI are log10 of the lowest and highest concentrations of a titration from 1e-10 to 1e-6.
"""

from __future__ import annotations

import math
import statistics

import pytest

from replicate_combine import (
    Filters,
    ReplicateValue,
    build_per_replicate_frame,
    collect_values,
    combine_replicates,
    merge_all,
    replicate_value,
)

LO = -10.0
HI = -6.0


def _row(kd: float | None, *, cls: str = "Good", out_of_range: bool = False, se: float | None = None, r2=None):
    return {"kd": kd, "affinityClass": cls, "kdOutOfRange": out_of_range, "kdLog10Se": se, "r2": r2}


def _value(label: str, log10_kd: float, *, se: float | None = None, state="measured", cls="Good"):
    return ReplicateValue(label, state, log10_kd, se, cls)


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
        m = combine_replicates([_value("A", -8.0), ReplicateValue("B", "failed", None, None, "Failed")])
        assert m.log10_kd == pytest.approx(-8.0)
        assert m.counts["failed"] == 1
        assert m.counts["measured"] == 1


class TestErrorModels:
    POINTS = [-8.0, -7.5, -8.6]

    def test_the_default_error_is_the_sample_sd_over_root_n(self):
        m = combine_replicates([_value(str(i), p) for i, p in enumerate(self.POINTS)])
        assert m.log10_kd_se == pytest.approx(statistics.stdev(self.POINTS) / math.sqrt(3))

    def test_the_phillips_error_is_the_population_sd_over_root_n_plus_one(self):
        m = combine_replicates([_value(str(i), p) for i, p in enumerate(self.POINTS)], error_model="phillips")
        assert m.log10_kd_se == pytest.approx(statistics.pstdev(self.POINTS) / math.sqrt(4))

    def test_the_two_error_models_report_the_same_value(self):
        values = [_value(str(i), p) for i, p in enumerate(self.POINTS)]
        assert combine_replicates(values).log10_kd == pytest.approx(
            combine_replicates(values, error_model="phillips").log10_kd
        )

    def test_one_replicate_keeps_its_own_fit_error_by_default(self):
        m = combine_replicates([_value("A", -8.0, se=0.31)])
        assert m.log10_kd_se == 0.31
        assert m.log10_kd_sd is None

    def test_one_replicate_with_no_fit_error_reports_none_and_keeps_its_own_class(self):
        m = combine_replicates([_value("A", -8.0, se=None, cls="Partial")])
        assert m.log10_kd_se is None
        assert m.classification.affinity_class == "Partial"

    def test_one_replicate_under_phillips_reports_exactly_zero(self):
        m = combine_replicates([_value("A", -8.0, se=0.31)], error_model="phillips")
        assert m.log10_kd_se == 0.0
        assert m.classification.affinity_class == "Good"

    def test_the_sd_column_is_the_sample_sd_whatever_the_error_model(self):
        values = [_value(str(i), p) for i, p in enumerate(self.POINTS)]
        expected = statistics.stdev(self.POINTS)
        assert combine_replicates(values).log10_kd_sd == pytest.approx(expected)
        assert combine_replicates(values, error_model="phillips").log10_kd_sd == pytest.approx(expected)

    def test_a_wide_spread_rates_partial(self):
        m = combine_replicates([_value("A", -9.0), _value("B", -6.0)])
        assert m.log10_kd_se > 0.5
        assert m.classification.affinity_class == "Partial"


class TestCensoring:
    def test_every_replicate_censored_low_reports_that_wall_with_no_error(self):
        m = combine_replicates([_value("A", LO, state="below"), _value("B", LO, state="below")])
        assert (m.log10_kd, m.log10_kd_se, m.out_of_range) == (LO, None, True)
        assert m.classification.affinity_class == "Partial"

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

    def test_a_majority_threshold_of_one_refuses_any_disagreement(self):
        m = combine_replicates(
            [_value("A", LO, state="below"), _value("B", LO, state="below"), _value("C", HI, state="above")],
            censored_majority=1.0,
        )
        assert m.classification.failure_reason == "conflicting_censoring"

    def test_one_measured_replicate_beside_opposite_walls_still_takes_the_mean(self):
        m = combine_replicates([_value("A", -8.0), _value("B", LO, state="below"), _value("C", HI, state="above")])
        assert m.log10_kd == pytest.approx((-8.0 + LO + HI) / 3)
        assert m.out_of_range is False


class TestFilters:
    def test_no_filter_is_active_by_default(self):
        assert Filters().active is False
        assert replicate_value("A", _row(1e-8, se=9.9, r2=0.01), log10_lo=LO, log10_hi=HI).state == "measured"

    def test_the_fit_error_cutoff_drops_an_imprecise_measured_replicate(self):
        filters = Filters(fit_error_cutoff=1.0)
        assert replicate_value("A", _row(1e-8, se=1.4), log10_lo=LO, log10_hi=HI, filters=filters).state == "filtered"
        assert replicate_value("A", _row(1e-8, se=0.9), log10_lo=LO, log10_hi=HI, filters=filters).state == "measured"

    def test_the_fit_error_cutoff_drops_a_measured_replicate_with_no_error(self):
        filters = Filters(fit_error_cutoff=1.0)
        assert replicate_value("A", _row(1e-8, se=None), log10_lo=LO, log10_hi=HI, filters=filters).state == "filtered"

    def test_a_censored_replicate_is_exempt_from_the_fit_error_cutoff(self):
        filters = Filters(fit_error_cutoff=0.01)
        v = replicate_value("A", _row(1e-14, out_of_range=True), log10_lo=LO, log10_hi=HI, filters=filters)
        assert v.state == "below"

    def test_the_r_squared_floor_drops_a_poor_fit(self):
        filters = Filters(r2_floor=0.8)
        assert (
            replicate_value("A", _row(1e-8, se=0.1, r2=0.5), log10_lo=LO, log10_hi=HI, filters=filters).state
            == "filtered"
        )
        assert (
            replicate_value("A", _row(1e-8, se=0.1, r2=0.9), log10_lo=LO, log10_hi=HI, filters=filters).state
            == "measured"
        )

    def test_a_filtered_replicate_is_left_out_of_the_mean_and_counted(self):
        values = [_value("A", -8.0), ReplicateValue("B", "filtered", -6.0, 2.0, "Partial")]
        m = combine_replicates(values)
        assert m.log10_kd == pytest.approx(-8.0)
        assert m.counts["filtered"] == 1

    def test_every_replicate_filtered_fails_with_a_named_reason(self):
        m = combine_replicates([ReplicateValue("A", "filtered", -8.0, 2.0, "Partial")])
        assert m.log10_kd is None
        assert m.classification.failure_reason == "all_replicates_filtered"

    def test_every_replicate_failed_names_no_filter(self):
        m = combine_replicates([ReplicateValue("A", "failed", None, None, "Failed")])
        assert m.classification.failure_reason is None


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
        assert frame.columns == ["clonotypeKey", "replicate", "kd", "kdLog10Se", "r2", "affinityClass", "kdState"]
        assert frame.height == 4
        assert frame.filter(kdState="failed").height == 2

    def test_the_counts_add_up_to_the_replicates_seen(self):
        merged = merge_all(self._collected())
        assert merged["v1"].counts == {"measured": 2, "below": 0, "above": 0, "failed": 0, "filtered": 0}
        assert merged["v2"].counts["failed"] == 2


class TestNonBinderAndRangeSide:
    def test_every_high_wall_replicate_a_non_binder_combines_as_a_non_binder(self):
        values = [_value(r, HI, state="above", cls="Non-binder") for r in ("A", "B")]
        merged = combine_replicates(values)
        assert merged.classification.affinity_class == "Non-binder"
        assert merged.range_side == "Above range"

    def test_one_high_wall_replicate_that_rose_keeps_the_variant_partial(self):
        values = [
            _value("A", HI, state="above", cls="Non-binder"),
            _value("B", HI, state="above", cls="Partial"),
        ]
        assert combine_replicates(values).classification.affinity_class == "Partial"

    def test_a_low_wall_bound_is_partial_below_range(self):
        merged = combine_replicates([_value(r, LO, state="below", cls="Partial") for r in ("A", "B")])
        assert (merged.classification.affinity_class, merged.range_side) == ("Partial", "Below range")

    def test_a_mean_of_measured_values_is_in_range(self):
        merged = combine_replicates([_value("A", -8.0, se=0.1), _value("B", -8.2, se=0.1)])
        assert merged.range_side == "In range"

    def test_no_value_has_no_side(self):
        merged = combine_replicates([_value("A", -8.0, state="failed", cls="Failed")])
        assert merged.range_side is None
