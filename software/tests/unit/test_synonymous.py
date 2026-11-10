"""The Starr mean over synonymous variants, its class, and the amino-acid parent delta."""

from __future__ import annotations

import math
import statistics

import polars as pl
import pytest

from constants import KD_SD_GOOD
from synonymous import (
    CLASS_ONE,
    CLASS_SEVERAL,
    CLASS_STOP,
    CLASS_SYNONYMOUS,
    CLASS_WILDTYPE,
    REASON_NO_UNMUTATED,
    REASON_SEVERAL_UNMUTATED,
    Member,
    SynonymousInputs,
    aa_class,
    combine,
    combine_synonymous,
    member_from_fit,
    nt_aa_variants,
    nt_class,
    parent_delta,
    resolve_nt_to_aa,
    resolve_parents,
    stop_variants,
)

LO, HI = -10.0, -6.0


def _fit_row(kd, cls="Good", out_of_range=False, se=None):
    return {"kd": kd, "affinityClass": cls, "kdOutOfRange": out_of_range, "kdLog10Se": se}


def _measured(value, se=None, cls="Good"):
    return Member("measured", value, se, cls)


class TestMemberFromFit:
    def test_a_measured_member_enters_at_its_own_log10_kd_with_its_se(self):
        c = member_from_fit(_fit_row(1e-8, se=0.1), LO, HI)
        assert c.state == "measured"
        assert c.log10_kd == pytest.approx(-8.0)
        assert c.log10_kd_se == 0.1

    def test_a_failed_member_is_left_out_even_when_it_reports_a_kd(self):
        assert member_from_fit(_fit_row(1e-8, cls="Failed"), LO, HI).state == "failed"

    def test_a_member_with_no_kd_is_failed(self):
        assert member_from_fit(_fit_row(None, cls="Partial"), LO, HI).state == "failed"

    def test_a_member_below_the_range_enters_at_the_lower_wall(self):
        c = member_from_fit(_fit_row(1e-13, out_of_range=True, se=0.2), LO, HI)
        assert (c.state, c.log10_kd, c.log10_kd_se) == ("below", LO, None)

    def test_a_member_above_the_range_enters_at_the_upper_wall(self):
        c = member_from_fit(_fit_row(1e-2, out_of_range=True), LO, HI)
        assert (c.state, c.log10_kd) == ("above", HI)

    def test_a_member_pinned_at_a_wall_keeps_its_side(self):
        assert member_from_fit(_fit_row(10.0**LO, out_of_range=True), LO, HI).state == "below"
        assert member_from_fit(_fit_row(10.0**HI, out_of_range=True), LO, HI).state == "above"


class TestCombine:
    def test_the_value_is_the_mean_of_log10_kd_and_the_se_is_sd_over_root_n(self):
        c = combine([_measured(-8.0), _measured(-7.0), _measured(-9.0)])
        assert c.log10_kd == pytest.approx(-8.0)
        assert c.log10_kd_se == pytest.approx(1.0 / math.sqrt(3))
        assert c.counts == {"measured": 3, "below": 0, "above": 0, "failed": 0}

    def test_the_fit_error_term_floors_the_se_when_the_members_agree_better_than_their_fits(self):
        c = combine([_measured(-8.0, se=0.3), _measured(-8.2, se=0.4)])
        # spread 0.1; fit-error term RMS(0.3, 0.4) / sqrt(2) = sqrt(0.125) / sqrt(2) = 0.25
        assert c.log10_kd_se == pytest.approx(0.25)
        assert c.classification.affinity_class == "Good"

    def test_the_spread_stands_when_it_exceeds_the_fit_error_term(self):
        c = combine([_measured(-7.0, se=0.1), _measured(-9.0, se=0.1)])
        assert c.log10_kd_se == pytest.approx(1.0)

    def test_a_member_without_a_fit_se_counts_at_the_typical_error_of_the_others(self):
        # RMS over the one fit SE there is, 0.6; n counts all three members in the mean.
        c = combine([_measured(-8.0, se=0.6), _measured(-8.0), Member("above", -8.0, None, "Partial")])
        assert c.log10_kd_se == pytest.approx(0.6 / math.sqrt(3))

    def test_a_member_at_a_wall_never_lowers_the_se(self):
        # Two measured members, then a third entering at the wall with no fit SE: the fit-error term
        # stays RMS / sqrt(n), so the added member cannot shrink it below what the fits support.
        two = combine([_measured(-6.1, se=0.6), _measured(-6.1, se=0.6)])
        three = combine([_measured(-6.1, se=0.6), _measured(-6.1, se=0.6), Member("above", -6.0, None, "Partial")])
        assert two.log10_kd_se == pytest.approx(0.6 / math.sqrt(2))
        assert three.log10_kd_se == pytest.approx(0.6 / math.sqrt(3))
        assert three.log10_kd_se > 0.6 * math.sqrt(2) / 3

    def test_no_member_with_a_fit_se_leaves_the_spread(self):
        c = combine([_measured(-8.0), _measured(-8.2)])
        assert c.log10_kd_se == pytest.approx(0.1)

    def test_with_no_fit_se_good_needs_more_than_half_the_members_rated_good(self):
        # The spread is the only evidence; tight agreement of members rated Partial is not Good.
        partial = combine([_measured(-8.0, cls="Partial"), _measured(-8.05, cls="Partial")])
        assert partial.classification.affinity_class == "Partial"
        half = combine([_measured(-8.0, cls="Good"), _measured(-8.05, cls="Partial")])
        assert half.classification.affinity_class == "Partial"
        most = combine([_measured(-8.0), _measured(-8.05), _measured(-8.02, cls="Partial")])
        assert most.classification.affinity_class == "Good"

    def test_the_good_majority_rule_does_not_apply_when_a_member_reports_a_fit_se(self):
        c = combine([_measured(-8.0, se=0.1, cls="Partial"), _measured(-8.05, cls="Partial")])
        assert c.classification.affinity_class == "Good"

    def test_the_floor_can_lift_a_good_spread_to_partial(self):
        c = combine([_measured(-8.0, se=1.2), _measured(-8.0, se=1.2)])
        assert c.log10_kd_se == pytest.approx(math.sqrt(2 * 1.2**2) / 2)
        assert c.classification.affinity_class == "Partial"

    def test_tight_agreement_rates_good_and_scatter_rates_partial(self):
        assert combine([_measured(-8.0), _measured(-8.1)]).classification.affinity_class == "Good"
        assert combine([_measured(-6.5), _measured(-9.5)]).classification.affinity_class == "Partial"

    def test_the_good_threshold_is_inclusive(self):
        # One synonymous variant carries its fit SE through unchanged, so the threshold is hit exactly.
        assert combine([_measured(-8.0, se=KD_SD_GOOD, cls="Partial")]).classification.affinity_class == "Good"

    def test_a_failed_member_is_counted_but_not_averaged(self):
        c = combine([_measured(-8.0), _measured(-8.2), Member("failed", None, None, "Failed")])
        assert c.log10_kd == pytest.approx(-8.1)
        assert c.counts["failed"] == 1

    def test_no_member_in_the_mean_is_failed_with_no_value(self):
        c = combine([Member("failed", None, None, "Failed")])
        assert c.classification.affinity_class == "Failed"
        assert c.classification.failure_reason == "all_variants_failed"
        assert (c.log10_kd, c.log10_kd_se, c.out_of_range) == (None, None, None)

    def test_members_censored_on_one_side_give_a_bound_with_no_se(self):
        c = combine([Member("below", LO, None, "Below range"), Member("below", LO, None, "Below range")])
        assert c.log10_kd == LO
        assert c.log10_kd_se is None
        assert c.out_of_range is True
        assert c.classification.affinity_class == "Below range"
        assert c.counts["below"] == 2

    def test_a_bound_is_rated_by_its_wall_under_the_log_fluorescence_fit(self):
        c = combine([Member("above", HI, None, "Above range"), Member("above", HI, None, "Above range")])
        assert c.classification.affinity_class == "Above range"

    def test_under_the_mean_bin_fit_a_bound_keeps_the_class_its_members_share(self):
        # The mean bin fit declares no wall class; its synonymous variants are rated by R².
        same = combine([Member("above", HI, None, "Good"), Member("above", HI, None, "Good")], wall_classes=False)
        assert same.classification.affinity_class == "Good"
        mixed = combine([Member("above", HI, None, "Good"), Member("above", HI, None, "Partial")], wall_classes=False)
        assert mixed.classification.affinity_class == "Partial"

    def test_members_at_opposite_walls_with_no_measured_member_fail(self):
        # Their mean would sit inside the range, a value no synonymous variant measured.
        c = combine([Member("below", LO, None, "Below range"), Member("above", HI, None, "Above range")])
        assert (c.log10_kd, c.log10_kd_se, c.out_of_range) == (None, None, None)
        assert c.classification.affinity_class == "Failed"
        assert c.classification.failure_reason == "conflicting_censoring"
        assert c.counts == {"measured": 0, "below": 1, "above": 1, "failed": 0}

    def test_a_majority_on_one_wall_sets_the_bound(self):
        members = [Member("below", LO, None, "Below range")] * 2 + [Member("above", HI, None, "Above range")]
        c = combine(members)
        assert (c.log10_kd, c.log10_kd_se, c.out_of_range) == (LO, None, True)
        assert c.classification.affinity_class == "Below range"

    def test_under_the_mean_bin_fit_a_majority_bound_reads_the_class_of_its_own_wall(self):
        members = [Member("above", HI, None, "Good")] * 2 + [Member("below", LO, None, "Partial")]
        c = combine(members, wall_classes=False)
        assert (c.log10_kd, c.classification.affinity_class) == (HI, "Good")

    def test_a_censored_member_beside_a_measured_one_enters_at_the_wall(self):
        c = combine([Member("above", HI, None, "Partial"), _measured(-8.0)])
        assert c.log10_kd == pytest.approx((HI - 8.0) / 2)
        assert c.out_of_range is False
        assert c.counts == {"measured": 1, "below": 0, "above": 1, "failed": 0}

    def test_one_member_with_a_fit_se_is_rated_on_that_se(self):
        assert combine([_measured(-8.0, se=0.2, cls="Partial")]).classification.affinity_class == "Good"
        c = combine([_measured(-8.0, se=0.9, cls="Good")])
        assert c.log10_kd_se == 0.9
        assert c.classification.affinity_class == "Partial"

    def test_one_member_with_no_se_keeps_its_own_class_and_a_null_se(self):
        c = combine([_measured(-8.0, se=None, cls="Good")])
        assert c.log10_kd_se is None
        assert c.classification.affinity_class == "Good"
        assert combine([_measured(-8.0, cls="Partial")]).classification.affinity_class == "Partial"


class TestSynonymousSd:
    def test_the_sd_is_the_synonymous_variants_sample_sd(self):
        c = combine([_measured(-8.0), _measured(-7.0), _measured(-9.0)])
        assert c.log10_kd_sd == pytest.approx(1.0)

    def test_the_sd_is_not_floored_by_the_fit_errors(self):
        c = combine([_measured(-8.0, se=0.3), _measured(-8.2, se=0.4)])
        assert c.log10_kd_se == pytest.approx(0.25)
        assert c.log10_kd_sd == pytest.approx(statistics.stdev([-8.0, -8.2]))

    def test_censored_members_enter_the_sd_at_their_wall(self):
        c = combine([_measured(-8.0), Member("above", HI, None, "Partial")])
        assert c.log10_kd_sd == pytest.approx(statistics.stdev([-8.0, HI]))

    def test_one_member_has_no_sd_even_with_a_fit_se(self):
        c = combine([_measured(-8.0, se=0.2), Member("failed", None, None, "Failed")])
        assert c.log10_kd_se == 0.2
        assert c.log10_kd_sd is None

    def test_a_bound_has_no_sd(self):
        c = combine([Member("below", LO, None, "Partial"), Member("below", LO, None, "Partial")])
        assert c.log10_kd_sd is None


class TestVariantClass:
    @pytest.mark.parametrize(
        ("mutations", "stop", "expected"),
        [
            (0, False, CLASS_WILDTYPE),
            (1, False, CLASS_ONE),
            (2, False, CLASS_SEVERAL),
            (7, False, CLASS_SEVERAL),
            (0, True, CLASS_STOP),
            (3, True, CLASS_STOP),
            (None, True, CLASS_STOP),
            (None, False, None),
        ],
    )
    def test_the_aa_class(self, mutations, stop, expected):
        assert aa_class(mutations, stop) == expected

    @pytest.mark.parametrize(
        ("aa_mutations", "nt_mutations", "stop", "expected"),
        [
            (0, 0, False, CLASS_WILDTYPE),
            (0, 2, False, CLASS_SYNONYMOUS),
            (1, 3, False, CLASS_ONE),
            (2, 2, False, CLASS_SEVERAL),
            (1, 1, True, CLASS_STOP),
            (0, None, False, None),
            (1, None, False, CLASS_ONE),
            (None, 0, False, None),
        ],
    )
    def test_the_nt_class_follows_its_aa_variant(self, aa_mutations, nt_mutations, stop, expected):
        assert nt_class(aa_mutations, nt_mutations, stop) == expected

    def test_only_false_marks_a_stop(self):
        productive = pl.DataFrame({"aaKey": ["a", "b", "c", "d"], "productive": ["false", "true", None, " False "]})
        assert stop_variants(productive) == {"a", "d"}


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
        failed = combine([Member("failed", None, None, "Failed")])
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
        frame, _, summary = combine_synonymous(
            self.PER_CLONOTYPE, SynonymousInputs(self.LINKS), log10_lo=LO, log10_hi=HI
        )
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
        # WT's two measured synonymous variants agree; M has one measured synonymous variant, so no spread.
        assert frame["kdLog10SynonymSd"].to_list() == [None, 0.0]

    def test_two_identical_members_have_the_fit_error_term_as_their_se(self):
        frame, _, _ = combine_synonymous(self.PER_CLONOTYPE, SynonymousInputs(self.LINKS), log10_lo=LO, log10_hi=HI)
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
        frame, _, summary = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        rows = {r["aaKey"]: r for r in frame.iter_rows(named=True)}
        assert (rows["WT"]["kdLog10ParentDelta"], rows["WT"]["kdLog10ParentDeltaSe"]) == (0.0, 0.0)
        assert rows["M"]["kdLog10ParentDelta"] == pytest.approx(-1.0)
        # M has one synonymous variant with a fit SE of 0.2; WT's two identical synonymous variants carry
        # their fit-error term.
        assert rows["M"]["kdLog10ParentDeltaSe"] == pytest.approx(math.hypot(0.2, math.sqrt(0.02) / 2))
        assert summary["parents"] == 1
        assert summary["unidentifiableParents"] == []
        assert summary["multiParentVariants"] == 0

    CLASS_INPUTS = {
        "aa_mutations": pl.DataFrame({"aaKey": ["WT", "M"], "mutationCount": ["0", "1"]}),
        "aa_parent_links": pl.DataFrame({"aaKey": ["WT", "M"], "parentId": ["P", "P"]}),
        "nt_mutations": pl.DataFrame({"ntKey": ["wt1", "wt2", "m1", "m2"], "mutationCount": ["0", "1", "1", "2"]}),
    }

    def test_each_grain_gets_its_class(self):
        inputs = SynonymousInputs(self.LINKS, **self.CLASS_INPUTS)
        frame, per_nt, _ = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        assert dict(zip(frame["aaKey"], frame["variantClass"])) == {"M": CLASS_ONE, "WT": CLASS_WILDTYPE}
        assert dict(zip(per_nt["ntKey"], per_nt["variantClass"])) == {
            "wt1": CLASS_WILDTYPE,
            "wt2": CLASS_SYNONYMOUS,
            "m1": CLASS_ONE,
            "m2": CLASS_ONE,
            "amb": None,
            "orphan": None,
        }

    def test_a_stop_variant_loses_its_values_but_keeps_its_counts(self):
        productive = pl.DataFrame({"aaKey": ["WT", "M"], "productive": ["true", "false"]})
        inputs = SynonymousInputs(self.LINKS, aa_productive=productive, **self.CLASS_INPUTS)
        frame, per_nt, summary = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        m = frame.filter(pl.col("aaKey") == "M").row(0, named=True)
        assert (m["kd"], m["kdLog10Se"], m["kdOutOfRange"]) == (None, None, None)
        assert (m["kdLog10ParentDelta"], m["kdLog10ParentDeltaSe"]) == (None, None)
        assert m["affinityClass"] == "Failed"
        assert m["fitFailureReason"] == "stop_codon"
        assert m["variantClass"] == CLASS_STOP
        assert (m["ntVariants"], m["ntMeasured"], m["ntFailed"]) == (2, 1, 1)
        assert summary["stopAaVariants"] == 1
        assert dict(zip(per_nt["ntKey"], per_nt["variantClass"]))["m1"] == CLASS_STOP
        wt = frame.filter(pl.col("aaKey") == "WT").row(0, named=True)
        assert wt["kd"] == pytest.approx(1e-8)

    def test_the_frame_names_why_an_amino_acid_variant_failed(self):
        per_clonotype = _per_clonotype(
            [
                ("lo", 10.0**LO, "Below range", True, None),
                ("hi", 10.0**HI, "Above range", True, None),
                ("bad", None, "Failed", None, None),
                ("ok", 1e-8, "Good", False, 0.1),
            ]
        )
        links = pl.DataFrame({"ntKey": ["lo", "hi", "bad", "ok"], "aaKey": ["CONFLICT", "CONFLICT", "NONE", "OK"]})
        frame, _, _ = combine_synonymous(per_clonotype, SynonymousInputs(links), log10_lo=LO, log10_hi=HI)
        assert dict(zip(frame["aaKey"], frame["fitFailureReason"])) == {
            "CONFLICT": "conflicting_censoring",
            "NONE": "all_variants_failed",
            "OK": None,
        }

    def test_without_the_productive_column_there_is_no_stop(self):
        inputs = SynonymousInputs(self.LINKS, **self.CLASS_INPUTS)
        frame, per_nt, summary = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        assert CLASS_STOP not in frame["variantClass"].to_list()
        assert CLASS_STOP not in per_nt["variantClass"].to_list()
        assert "stopAaVariants" not in summary
        assert frame.filter(pl.col("aaKey") == "M")["kd"].item() == pytest.approx(1e-9)

    def test_without_the_aa_mutation_count_there_is_no_class_on_either_grain(self):
        productive = pl.DataFrame({"aaKey": ["M"], "productive": ["false"]})
        inputs = SynonymousInputs(self.LINKS, aa_productive=productive, nt_mutations=self.CLASS_INPUTS["nt_mutations"])
        frame, per_nt, summary = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        assert "variantClass" not in frame.columns
        assert per_nt.columns == ["ntKey"]
        assert frame.filter(pl.col("aaKey") == "M")["kd"].item() is None
        assert summary["stopAaVariants"] == 1

    def test_the_aa_class_needs_no_parent_link(self):
        inputs = SynonymousInputs(self.LINKS, aa_mutations=self.CLASS_INPUTS["aa_mutations"])
        frame, per_nt, summary = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        assert "variantClass" in frame.columns
        assert "kdLog10ParentDelta" not in frame.columns
        assert per_nt is None
        assert summary["parentMode"] is False

    NT_PARENT_LINKS = pl.DataFrame({"ntKey": ["wt1", "wt2", "m1", "m2"], "parentId": ["P", "P", "P", "P"]})

    def test_the_nt_parent_delta_reads_each_nt_variant_against_the_unmutated_one(self):
        per_clonotype = _per_clonotype(
            [
                ("wt1", 1e-8, "Good", False, 0.1),
                ("wt2", 1e-7, "Good", False, 0.2),
                ("m1", 1e-9, "Good", False, None),
                ("m2", None, "Failed", None, None),
            ]
        )
        inputs = SynonymousInputs(
            self.LINKS,
            nt_mutations=self.CLASS_INPUTS["nt_mutations"],
            nt_parent_links=self.NT_PARENT_LINKS,
        )
        _, per_nt, summary = combine_synonymous(per_clonotype, inputs, log10_lo=LO, log10_hi=HI)
        rows = {r["ntKey"]: r for r in per_nt.iter_rows(named=True)}
        assert (rows["wt1"]["kdLog10ParentDelta"], rows["wt1"]["kdLog10ParentDeltaSe"]) == (0.0, 0.0)
        assert rows["wt2"]["kdLog10ParentDelta"] == pytest.approx(1.0)
        assert rows["wt2"]["kdLog10ParentDeltaSe"] == pytest.approx(math.hypot(0.1, 0.2))
        assert rows["m1"]["kdLog10ParentDelta"] == pytest.approx(-1.0)
        assert rows["m1"]["kdLog10ParentDeltaSe"] is None
        assert (rows["m2"]["kdLog10ParentDelta"], rows["m2"]["kdLog10ParentDeltaSe"]) == (None, None)
        assert "variantClass" not in per_nt.columns
        assert summary["ntParentMode"] is True
        assert (summary["ntParents"], summary["ntUnidentifiableParents"], summary["ntMultiParentVariants"]) == (
            1,
            [],
            0,
        )

    def test_a_censored_nt_variant_reads_from_its_wall(self):
        per_clonotype = _per_clonotype([("wt1", 1e-8, "Good", False, 0.1), ("m1", 1e-13, "Partial", True, None)])
        inputs = SynonymousInputs(
            self.LINKS,
            nt_mutations=self.CLASS_INPUTS["nt_mutations"],
            nt_parent_links=self.NT_PARENT_LINKS,
        )
        _, per_nt, _ = combine_synonymous(per_clonotype, inputs, log10_lo=LO, log10_hi=HI)
        m1 = per_nt.filter(pl.col("ntKey") == "m1").row(0, named=True)
        assert m1["kdLog10ParentDelta"] == pytest.approx(LO - -8.0)
        assert m1["kdLog10ParentDeltaSe"] is None

    def test_an_nt_parent_with_none_or_several_unmutated_nt_variants_is_not_identifiable(self):
        per_clonotype = _per_clonotype([(nt, 1e-8, "Good", False, 0.1) for nt in ("wt1", "wt2", "m1", "m2")])
        links = pl.DataFrame({"ntKey": ["wt1", "wt2", "m1", "m2"], "parentId": ["P", "P", "Q", "Q"]})
        mutations = pl.DataFrame({"ntKey": ["wt1", "wt2", "m1", "m2"], "mutationCount": ["0", "0", "1", "2"]})
        inputs = SynonymousInputs(self.LINKS, nt_mutations=mutations, nt_parent_links=links)
        _, per_nt, summary = combine_synonymous(per_clonotype, inputs, log10_lo=LO, log10_hi=HI)
        assert per_nt["kdLog10ParentDelta"].to_list() == [None, None, None, None]
        assert summary["ntUnidentifiableParents"] == [
            {"parentId": "P", "reason": REASON_SEVERAL_UNMUTATED},
            {"parentId": "Q", "reason": REASON_NO_UNMUTATED},
        ]

    def test_a_stop_nt_variant_keeps_its_nt_parent_delta(self):
        productive = pl.DataFrame({"aaKey": ["WT", "M"], "productive": ["true", "false"]})
        inputs = SynonymousInputs(
            self.LINKS,
            aa_productive=productive,
            nt_mutations=self.CLASS_INPUTS["nt_mutations"],
            nt_parent_links=self.NT_PARENT_LINKS,
        )
        _, per_nt, _ = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        assert per_nt.filter(pl.col("ntKey") == "m1")["kdLog10ParentDelta"].item() == pytest.approx(-1.0)

    def test_nt_parent_links_without_nt_mutation_counts_are_refused(self):
        with pytest.raises(ValueError, match="nucleotide mutation counts"):
            SynonymousInputs(self.LINKS, nt_parent_links=self.NT_PARENT_LINKS)

    def test_parent_links_without_mutation_counts_are_refused(self):
        with pytest.raises(ValueError, match="mutation counts"):
            SynonymousInputs(self.LINKS, aa_parent_links=pl.DataFrame({"aaKey": [], "parentId": []}))


class TestPartlyCensoredMembers:
    def test_members_at_a_wall_cap_the_class_at_partial(self):
        c = combine([Member("above", -5.0, None, "Good"), Member("above", -5.0, None, "Good"), _measured(-5.9, se=0.2)])
        assert c.log10_kd_se is not None and c.log10_kd_se <= 0.5
        assert c.classification.affinity_class == "Partial"
        assert c.out_of_range is False

    def test_measured_members_alone_still_rate_good(self):
        assert combine([_measured(-8.0, se=0.1), _measured(-8.1, se=0.1)]).classification.affinity_class == "Good"

    def test_a_measured_majority_rates_good_with_one_member_at_a_wall(self):
        c = combine([_measured(-9.8, se=0.1), _measured(-9.9, se=0.1), Member("below", -10.0, None, "Good")])
        assert c.classification.affinity_class == "Good"

    def test_half_measured_members_are_not_a_majority(self):
        c = combine(
            [
                _measured(-9.8, se=0.1),
                _measured(-9.9, se=0.1),
                Member("below", -10.0, None, "Good"),
                Member("below", -10.0, None, "Good"),
            ]
        )
        assert c.log10_kd_se is not None and c.log10_kd_se <= 0.5
        assert c.classification.affinity_class == "Partial"


class TestNtAaVariants:
    # A and B have 2 synonymous nt variants each; C has 1; n4 maps to both A and B.
    LINKS = pl.DataFrame(
        {
            "ntKey": ["n1", "n2", "n3", "n5", "n6", "n4", "n4"],
            "aaKey": ["A", "A", "B", "B", "C", "A", "B"],
        }
    )

    @staticmethod
    def _variants(frame: pl.DataFrame) -> dict[str, str | None]:
        return dict(zip(frame["clonotypeKey"], frame["aaVariant"]))

    def test_each_nt_variant_names_its_amino_acid_variant_by_label(self):
        labels = pl.DataFrame({"aaKey": ["A", "B"], "label": ["V-ALPHA", "V-BETA"]})
        assert self._variants(nt_aa_variants(self.LINKS, labels)) == {
            "n1": "V-ALPHA",
            "n2": "V-ALPHA",
            "n3": "V-BETA",
            "n5": "V-BETA",
            "n6": None,
        }

    @pytest.mark.parametrize("order", [["V-ZETA", "V-ALPHA"], ["V-ALPHA", "V-ZETA"]])
    def test_an_amino_acid_variant_with_two_labels_takes_the_smallest(self, order):
        labels = pl.DataFrame({"aaKey": ["A", "A", "B"], "label": [*order, "V-BETA"]})
        for _ in range(5):
            assert self._variants(nt_aa_variants(self.LINKS, labels))["n1"] == "V-ALPHA"

    def test_a_variant_with_no_label_keeps_its_key(self):
        frame = nt_aa_variants(self.LINKS, pl.DataFrame({"aaKey": ["A"], "label": ["V-ALPHA"]}))
        assert self._variants(frame)["n3"] == "B"

    def test_every_single_mapped_nt_variant_carries_its_amino_acid_key(self):
        frame = nt_aa_variants(self.LINKS, None)
        assert dict(zip(frame["clonotypeKey"], frame["aaKey"])) == {
            "n1": "A",
            "n2": "A",
            "n3": "B",
            "n5": "B",
            "n6": "C",
        }

    def test_an_nt_variant_with_several_amino_acid_variants_is_left_out(self):
        assert "n4" not in nt_aa_variants(self.LINKS, None)["clonotypeKey"].to_list()

    def test_an_amino_acid_variant_with_one_nt_variant_gets_no_label(self):
        assert self._variants(nt_aa_variants(self.LINKS, None))["n6"] is None

    def test_only_nt_variants_with_a_curve_count_towards_the_two(self):
        # B keeps one curve, so neither of its nt variants is listed; A keeps both.
        variants = self._variants(nt_aa_variants(self.LINKS, None, curved=["n1", "n2", "n3"]))
        assert variants["n1"] == "A" and variants["n2"] == "A"
        assert variants["n3"] is None and variants["n5"] is None

    def test_an_nt_variant_without_a_curve_gets_no_label(self):
        variants = self._variants(nt_aa_variants(self.LINKS, None, curved=["n1", "n2", "n3", "n5", "n6"]))
        assert variants["n6"] is None
