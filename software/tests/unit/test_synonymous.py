"""The Starr mean over synonymous copies, its class, and the amino-acid parent delta."""

from __future__ import annotations

import math

import polars as pl
import pytest

from constants import KD_SD_GOOD
from synonymous import (
    REASON_NO_UNMUTATED,
    REASON_SEVERAL_UNMUTATED,
    Copy,
    SynonymousInputs,
    combine,
    combine_synonymous,
    copy_from_fit,
    parent_delta,
    resolve_nt_to_aa,
    resolve_parents,
)

LO, HI = -10.0, -6.0


def _fit_row(kd, cls="Good", out_of_range=False, se=None):
    return {"kd": kd, "affinityClass": cls, "kdOutOfRange": out_of_range, "kdLog10Se": se}


def _measured(value, se=None, cls="Good"):
    return Copy("measured", value, se, cls)


class TestCopyFromFit:
    def test_a_measured_copy_enters_at_its_own_log10_kd_with_its_se(self):
        c = copy_from_fit(_fit_row(1e-8, se=0.1), LO, HI)
        assert c.state == "measured"
        assert c.log10_kd == pytest.approx(-8.0)
        assert c.log10_kd_se == 0.1

    def test_a_failed_copy_is_left_out_even_when_it_reports_a_kd(self):
        assert copy_from_fit(_fit_row(1e-8, cls="Failed"), LO, HI).state == "failed"

    def test_a_copy_with_no_kd_is_failed(self):
        assert copy_from_fit(_fit_row(None, cls="Partial"), LO, HI).state == "failed"

    def test_a_copy_below_the_range_enters_at_the_lower_wall(self):
        c = copy_from_fit(_fit_row(1e-13, out_of_range=True, se=0.2), LO, HI)
        assert (c.state, c.log10_kd, c.log10_kd_se) == ("below", LO, None)

    def test_a_copy_above_the_range_enters_at_the_upper_wall(self):
        c = copy_from_fit(_fit_row(1e-2, out_of_range=True), LO, HI)
        assert (c.state, c.log10_kd) == ("above", HI)

    def test_a_copy_pinned_at_a_wall_keeps_its_side(self):
        assert copy_from_fit(_fit_row(10.0**LO, out_of_range=True), LO, HI).state == "below"
        assert copy_from_fit(_fit_row(10.0**HI, out_of_range=True), LO, HI).state == "above"


class TestCombine:
    def test_the_value_is_the_mean_of_log10_kd_and_the_se_is_sd_over_root_n(self):
        c = combine([_measured(-8.0), _measured(-7.0), _measured(-9.0)])
        assert c.log10_kd == pytest.approx(-8.0)
        assert c.log10_kd_se == pytest.approx(1.0 / math.sqrt(3))
        assert c.counts == {"measured": 3, "below": 0, "above": 0, "failed": 0}

    def test_the_fit_error_term_floors_the_se_when_the_copies_agree_better_than_their_fits(self):
        c = combine([_measured(-8.0, se=0.3), _measured(-8.2, se=0.4)])
        # spread 0.1; fit-error term sqrt(0.3^2 + 0.4^2) / 2 = 0.25
        assert c.log10_kd_se == pytest.approx(0.25)
        assert c.classification.affinity_class == "Good"

    def test_the_spread_stands_when_it_exceeds_the_fit_error_term(self):
        c = combine([_measured(-7.0, se=0.1), _measured(-9.0, se=0.1)])
        assert c.log10_kd_se == pytest.approx(1.0)

    def test_a_copy_without_a_fit_se_is_left_out_of_the_sum_but_not_the_divisor(self):
        c = combine([_measured(-8.0, se=0.6), _measured(-8.0), Copy("above", -8.0, None, "Partial")])
        assert c.log10_kd_se == pytest.approx(0.6 / 3)

    def test_no_copy_with_a_fit_se_leaves_the_spread(self):
        c = combine([_measured(-8.0), _measured(-8.2)])
        assert c.log10_kd_se == pytest.approx(0.1)

    def test_the_floor_can_lift_a_good_spread_to_partial(self):
        c = combine([_measured(-8.0, se=1.2), _measured(-8.0, se=1.2)])
        assert c.log10_kd_se == pytest.approx(math.sqrt(2 * 1.2**2) / 2)
        assert c.classification.affinity_class == "Partial"

    def test_tight_agreement_rates_good_and_scatter_rates_partial(self):
        assert combine([_measured(-8.0), _measured(-8.1)]).classification.affinity_class == "Good"
        assert combine([_measured(-6.5), _measured(-9.5)]).classification.affinity_class == "Partial"

    def test_the_good_threshold_is_inclusive(self):
        # One copy carries its fit SE through unchanged, so the threshold is hit exactly.
        assert combine([_measured(-8.0, se=KD_SD_GOOD, cls="Partial")]).classification.affinity_class == "Good"

    def test_a_failed_copy_is_counted_but_not_averaged(self):
        c = combine([_measured(-8.0), _measured(-8.2), Copy("failed", None, None, "Failed")])
        assert c.log10_kd == pytest.approx(-8.1)
        assert c.counts["failed"] == 1

    def test_no_copy_in_the_mean_is_failed_with_no_value(self):
        c = combine([Copy("failed", None, None, "Failed")])
        assert c.classification.affinity_class == "Failed"
        assert (c.log10_kd, c.log10_kd_se, c.out_of_range) == (None, None, None)

    def test_copies_censored_on_one_side_give_a_bound_with_no_se(self):
        c = combine([Copy("below", LO, None, "Partial"), Copy("below", LO, None, "Partial")])
        assert c.log10_kd == LO
        assert c.log10_kd_se is None
        assert c.out_of_range is True
        assert c.classification.affinity_class == "Partial"
        assert c.counts["below"] == 2

    def test_a_censored_copy_beside_a_measured_one_enters_at_the_wall(self):
        c = combine([Copy("above", HI, None, "Partial"), _measured(-8.0)])
        assert c.log10_kd == pytest.approx((HI - 8.0) / 2)
        assert c.out_of_range is False
        assert c.counts == {"measured": 1, "below": 0, "above": 1, "failed": 0}

    def test_one_copy_with_a_fit_se_is_rated_on_that_se(self):
        assert combine([_measured(-8.0, se=0.2, cls="Partial")]).classification.affinity_class == "Good"
        c = combine([_measured(-8.0, se=0.9, cls="Good")])
        assert c.log10_kd_se == 0.9
        assert c.classification.affinity_class == "Partial"

    def test_one_copy_with_no_se_keeps_its_own_class_and_a_null_se(self):
        c = combine([_measured(-8.0, se=None, cls="Good")])
        assert c.log10_kd_se is None
        assert c.classification.affinity_class == "Good"
        assert combine([_measured(-8.0, cls="Partial")]).classification.affinity_class == "Partial"


class TestResolveNtToAa:
    def test_a_variant_with_one_aa_maps_and_the_rest_are_counted(self):
        links = pl.DataFrame({"ntKey": ["n1", "n2", "n3", "n3"], "aaKey": ["a1", "a1", "a1", "a2"]})
        mapping, unmapped, ambiguous = resolve_nt_to_aa(links, ["n1", "n2", "n3", "n4"])
        assert mapping == {"n1": "a1", "n2": "a1"}
        assert (unmapped, ambiguous) == (1, 1)

    def test_a_repeated_link_is_one_link(self):
        links = pl.DataFrame({"ntKey": ["n1", "n1"], "aaKey": ["a1", "a1"]})
        assert resolve_nt_to_aa(links, ["n1"]) == ({"n1": "a1"}, 0, 0)


class TestResolveParents:
    def test_the_parent_is_its_single_unmutated_variant(self):
        links = pl.DataFrame({"aaKey": ["wt", "m1", "m2"], "parentId": ["P", "P", "P"]})
        mutations = pl.DataFrame({"aaKey": ["wt", "m1", "m2"], "mutationCount": ["0", "1", "2"]})
        assignment, own, unidentifiable, several = resolve_parents(links, mutations)
        assert assignment == {"wt": "P", "m1": "P", "m2": "P"}
        assert own == {"P": "wt"}
        assert (unidentifiable, several) == ([], 0)

    def test_a_parent_with_none_or_several_unmutated_variants_is_not_identifiable(self):
        links = pl.DataFrame({"aaKey": ["a", "b", "c", "d"], "parentId": ["P", "P", "Q", "Q"]})
        mutations = pl.DataFrame({"aaKey": ["a", "b", "c", "d"], "mutationCount": ["1", "1", "0", "0"]})
        _, own, unidentifiable, _ = resolve_parents(links, mutations)
        assert own == {}
        assert unidentifiable == [
            {"parentId": "P", "reason": REASON_NO_UNMUTATED},
            {"parentId": "Q", "reason": REASON_SEVERAL_UNMUTATED},
        ]

    def test_a_variant_with_several_parents_has_none(self):
        links = pl.DataFrame({"aaKey": ["a", "a"], "parentId": ["P", "Q"]})
        mutations = pl.DataFrame({"aaKey": ["a"], "mutationCount": ["3"]})
        assignment, _, _, several = resolve_parents(links, mutations)
        assert assignment == {}
        assert several == 1


class TestParentDelta:
    def test_the_delta_is_the_difference_with_errors_in_quadrature(self):
        variant = combine([_measured(-9.0), _measured(-9.2)])
        parent = combine([_measured(-8.0), _measured(-8.6)])
        delta, se = parent_delta(variant, parent, is_parent=False)
        assert delta == pytest.approx(-9.1 - -8.3)
        assert se == pytest.approx(math.hypot(variant.log10_kd_se, parent.log10_kd_se))

    def test_the_parent_against_itself_is_zero_with_zero_error(self):
        parent = combine([_measured(-8.0), _measured(-8.6)])
        assert parent_delta(parent, parent, is_parent=True) == (0.0, 0.0)

    def test_a_missing_se_on_either_side_leaves_the_error_null(self):
        variant = combine([_measured(-9.0)])
        parent = combine([_measured(-8.0), _measured(-8.6)])
        delta, se = parent_delta(variant, parent, is_parent=False)
        assert delta == pytest.approx(-0.7)
        assert se is None

    def test_no_kd_on_either_side_gives_no_delta(self):
        failed = combine([Copy("failed", None, None, "Failed")])
        measured = combine([_measured(-8.0)])
        assert parent_delta(failed, measured, is_parent=False) == (None, None)
        assert parent_delta(measured, failed, is_parent=False) == (None, None)
        assert parent_delta(measured, None, is_parent=False) == (None, None)


def _per_clonotype(rows: list[tuple]) -> pl.DataFrame:
    return pl.DataFrame(
        rows,
        schema={
            "clonotypeKey": pl.Utf8,
            "kd": pl.Float64,
            "affinityClass": pl.Utf8,
            "kdOutOfRange": pl.Boolean,
            "kdLog10Se": pl.Float64,
        },
        orient="row",
    )


class TestCombineSynonymous:
    PER_CLONOTYPE = _per_clonotype(
        [
            ("wt1", 1e-8, "Good", False, 0.1),
            ("wt2", 1e-8, "Good", False, 0.1),
            ("m1", 1e-9, "Good", False, 0.2),
            ("m2", None, "Failed", None, None),
            ("amb", 1e-8, "Good", False, 0.1),
            ("orphan", 1e-8, "Good", False, 0.1),
        ]
    )
    LINKS = pl.DataFrame(
        {"ntKey": ["wt1", "wt2", "m1", "m2", "amb", "amb"], "aaKey": ["WT", "WT", "M", "M", "WT", "M"]}
    )

    def test_one_row_per_aa_variant_with_counts(self):
        frame, summary = combine_synonymous(self.PER_CLONOTYPE, SynonymousInputs(self.LINKS), log10_lo=LO, log10_hi=HI)
        assert frame["aaKey"].to_list() == ["M", "WT"]
        assert frame["ntVariants"].to_list() == [2, 2]
        assert frame["ntMeasured"].to_list() == [1, 2]
        assert frame["ntFailed"].to_list() == [1, 0]
        assert "kdLog10ParentDelta" not in frame.columns
        assert summary == {
            "ntVariants": 6,
            "aaVariants": 2,
            "unmappedNtVariants": 1,
            "ambiguousNtVariants": 1,
            "parentMode": False,
        }

    def test_two_identical_copies_have_the_fit_error_term_as_their_se(self):
        frame, _ = combine_synonymous(self.PER_CLONOTYPE, SynonymousInputs(self.LINKS), log10_lo=LO, log10_hi=HI)
        wt = frame.filter(pl.col("aaKey") == "WT").row(0, named=True)
        assert wt["kd"] == pytest.approx(1e-8)
        # zero scatter; the fit-error term is sqrt(0.1^2 + 0.1^2) / 2
        assert wt["kdLog10Se"] == pytest.approx(math.sqrt(0.02) / 2)
        assert wt["affinityClass"] == "Good"

    def test_the_parent_delta_is_added_when_the_parent_columns_come(self):
        inputs = SynonymousInputs(
            self.LINKS,
            aa_mutations=pl.DataFrame({"aaKey": ["WT", "M"], "mutationCount": ["0", "1"]}),
            aa_parent_links=pl.DataFrame({"aaKey": ["WT", "M"], "parentId": ["P", "P"]}),
        )
        frame, summary = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        rows = {r["aaKey"]: r for r in frame.iter_rows(named=True)}
        assert (rows["WT"]["kdLog10ParentDelta"], rows["WT"]["kdLog10ParentDeltaSe"]) == (0.0, 0.0)
        assert rows["M"]["kdLog10ParentDelta"] == pytest.approx(-1.0)
        # M has one copy with a fit SE of 0.2; WT's two identical copies carry their fit-error term.
        assert rows["M"]["kdLog10ParentDeltaSe"] == pytest.approx(math.hypot(0.2, math.sqrt(0.02) / 2))
        assert summary["parents"] == 1
        assert summary["unidentifiableParents"] == []
        assert summary["multiParentVariants"] == 0

    def test_half_the_parent_columns_is_refused(self):
        with pytest.raises(ValueError, match="together"):
            SynonymousInputs(self.LINKS, aa_mutations=pl.DataFrame({"aaKey": [], "mutationCount": []}))
