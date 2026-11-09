"""The Starr mean over synonymous copies, its class, and the amino-acid parent delta."""

from __future__ import annotations

import math
import random
import statistics

import polars as pl
import pytest

from classify import Classification
from constants import KD_SD_GOOD, SYNONYM_OUTLIER_Z
from replicate_combine import Merged
from synonymous import (
    CLASS_ONE,
    CLASS_SEVERAL,
    CLASS_STOP,
    CLASS_SYNONYMOUS,
    CLASS_WILDTYPE,
    OUTLIER_SCHEMA,
    REASON_NO_UNMUTATED,
    REASON_SEVERAL_UNMUTATED,
    Copy,
    OutlierTest,
    SynonymousInputs,
    aa_class,
    combine,
    combine_synonymous,
    copy_from_fit,
    dataset_sigma,
    nt_aa_variants,
    nt_class,
    outlier_groups,
    outlier_note,
    outlier_tests,
    parent_delta,
    resolve_nt_to_aa,
    resolve_parents,
    stop_variants,
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
        c = combine([Copy("below", LO, None, "Below range"), Copy("below", LO, None, "Below range")])
        assert c.log10_kd == LO
        assert c.log10_kd_se is None
        assert c.out_of_range is True
        assert c.classification.affinity_class == "Below range"
        assert c.counts["below"] == 2

    def test_a_bound_is_rated_by_its_wall_under_the_log_fluorescence_fit(self):
        c = combine([Copy("above", HI, None, "Above range"), Copy("above", HI, None, "Above range")])
        assert c.classification.affinity_class == "Above range"

    def test_under_the_mean_bin_fit_a_bound_keeps_the_class_its_copies_share(self):
        # The mean bin fit declares no wall class; its copies are rated by R².
        same = combine([Copy("above", HI, None, "Good"), Copy("above", HI, None, "Good")], wall_classes=False)
        assert same.classification.affinity_class == "Good"
        mixed = combine([Copy("above", HI, None, "Good"), Copy("above", HI, None, "Partial")], wall_classes=False)
        assert mixed.classification.affinity_class == "Partial"

    def test_copies_at_opposite_walls_with_no_measured_copy_fail(self):
        # Their mean would sit inside the range, a value no copy measured.
        c = combine([Copy("below", LO, None, "Below range"), Copy("above", HI, None, "Above range")])
        assert (c.log10_kd, c.log10_kd_se, c.out_of_range) == (None, None, None)
        assert c.classification.affinity_class == "Failed"
        assert c.classification.failure_reason == "conflicting_censoring"
        assert c.counts == {"measured": 0, "below": 1, "above": 1, "failed": 0}

    def test_a_majority_on_one_wall_sets_the_bound(self):
        copies = [Copy("below", LO, None, "Below range")] * 2 + [Copy("above", HI, None, "Above range")]
        c = combine(copies)
        assert (c.log10_kd, c.log10_kd_se, c.out_of_range) == (LO, None, True)
        assert c.classification.affinity_class == "Below range"

    def test_under_the_mean_bin_fit_a_majority_bound_reads_the_class_of_its_own_wall(self):
        copies = [Copy("above", HI, None, "Good")] * 2 + [Copy("below", LO, None, "Partial")]
        c = combine(copies, wall_classes=False)
        assert (c.log10_kd, c.classification.affinity_class) == (HI, "Good")

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


class TestSynonymSd:
    def test_the_sd_is_the_synonyms_sample_sd(self):
        c = combine([_measured(-8.0), _measured(-7.0), _measured(-9.0)])
        assert c.log10_kd_sd == pytest.approx(1.0)

    def test_the_sd_is_not_floored_by_the_fit_errors(self):
        c = combine([_measured(-8.0, se=0.3), _measured(-8.2, se=0.4)])
        assert c.log10_kd_se == pytest.approx(0.25)
        assert c.log10_kd_sd == pytest.approx(statistics.stdev([-8.0, -8.2]))

    def test_censored_copies_enter_the_sd_at_their_wall(self):
        c = combine([_measured(-8.0), Copy("above", HI, None, "Partial")])
        assert c.log10_kd_sd == pytest.approx(statistics.stdev([-8.0, HI]))

    def test_one_copy_has_no_sd_even_with_a_fit_se(self):
        c = combine([_measured(-8.0, se=0.2), Copy("failed", None, None, "Failed")])
        assert c.log10_kd_se == 0.2
        assert c.log10_kd_sd is None

    def test_a_bound_has_no_sd(self):
        c = combine([Copy("below", LO, None, "Partial"), Copy("below", LO, None, "Partial")])
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
            "synonymSpread": {"variants": 1, "median": 0.0},
            "ntParentMode": False,
            # WT's two copies are tested and agree; M has one copy in its mean.
            "synonymOutlier": {
                "zThreshold": 2.5,
                "datasetSigma": 0.0,
                "testedAaVariants": 1,
                "testedNtVariants": 2,
                "flaggedNtVariants": 0,
                "flaggedAaVariants": 0,
            },
        }
        # WT's two measured copies agree; M has one measured copy, so no spread.
        assert frame["kdLog10SynonymSd"].to_list() == [None, 0.0]

    def test_two_identical_copies_have_the_fit_error_term_as_their_se(self):
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
        # M has one copy with a fit SE of 0.2; WT's two identical copies carry their fit-error term.
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
        assert m["variantClass"] == CLASS_STOP
        assert (m["ntVariants"], m["ntMeasured"], m["ntFailed"]) == (2, 1, 1)
        assert summary["stopAaVariants"] == 1
        assert dict(zip(per_nt["ntKey"], per_nt["variantClass"]))["m1"] == CLASS_STOP
        wt = frame.filter(pl.col("aaKey") == "WT").row(0, named=True)
        assert wt["kd"] == pytest.approx(1e-8)

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
        assert per_nt.columns == ["ntKey", *OUTLIER_SCHEMA]
        assert frame.filter(pl.col("aaKey") == "M")["kd"].item() is None
        assert summary["stopAaVariants"] == 1

    def test_the_aa_class_needs_no_parent_link(self):
        inputs = SynonymousInputs(self.LINKS, aa_mutations=self.CLASS_INPUTS["aa_mutations"])
        frame, per_nt, summary = combine_synonymous(self.PER_CLONOTYPE, inputs, log10_lo=LO, log10_hi=HI)
        assert "variantClass" in frame.columns
        assert "kdLog10ParentDelta" not in frame.columns
        assert per_nt.columns == ["ntKey", *OUTLIER_SCHEMA]
        assert summary["parentMode"] is False

    NT_PARENT_LINKS = pl.DataFrame({"ntKey": ["wt1", "wt2", "m1", "m2"], "parentId": ["P", "P", "P", "P"]})

    def test_the_nt_parent_delta_reads_each_copy_against_the_unmutated_copy(self):
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

    def test_a_censored_copy_reads_from_its_wall(self):
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

    def test_an_nt_parent_with_none_or_several_unmutated_copies_is_not_identifiable(self):
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

    def test_a_stop_copy_keeps_its_nt_parent_delta(self):
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


class TestPartlyCensoredCopies:
    def test_copies_at_a_wall_cap_the_class_at_partial(self):
        c = combine([Copy("above", -5.0, None, "Good"), Copy("above", -5.0, None, "Good"), _measured(-5.9, se=0.2)])
        assert c.log10_kd_se is not None and c.log10_kd_se <= 0.5
        assert c.classification.affinity_class == "Partial"
        assert c.out_of_range is False

    def test_measured_copies_alone_still_rate_good(self):
        assert combine([_measured(-8.0, se=0.1), _measured(-8.1, se=0.1)]).classification.affinity_class == "Good"

    def test_a_measured_majority_rates_good_with_one_copy_at_a_wall(self):
        c = combine([_measured(-9.8, se=0.1), _measured(-9.9, se=0.1), Copy("below", -10.0, None, "Good")])
        assert c.classification.affinity_class == "Good"

    def test_half_measured_copies_are_not_a_majority(self):
        c = combine(
            [
                _measured(-9.8, se=0.1),
                _measured(-9.9, se=0.1),
                Copy("below", -10.0, None, "Good"),
                Copy("below", -10.0, None, "Good"),
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


def _outlier_run(groups: dict[str, list[tuple[float, float | None]]]):
    """combine_synonymous over measured copies given as (log10 K_D, fit SE); copy i of `aa` is `aa`+i."""
    rows, links = [], []
    for aa, copies in groups.items():
        for i, (log10_kd, se) in enumerate(copies, start=1):
            rows.append((f"{aa}{i}", 10.0**log10_kd, "Good", False, se))
            links.append((f"{aa}{i}", aa))
    inputs = SynonymousInputs(pl.DataFrame(links, schema=["ntKey", "aaKey"], orient="row"))
    frame, per_nt, summary = combine_synonymous(_per_clonotype(rows), inputs, log10_lo=LO, log10_hi=HI)
    return frame, {r["ntKey"]: r for r in per_nt.iter_rows(named=True)}, summary


def _note(frame: pl.DataFrame, aa: str) -> str | None:
    return frame.filter(pl.col("aaKey") == aa)["synonymOutliers"].item()


class TestSynonymOutlier:
    # Hand-checked: delta against the mean of the other two; the others' SE is max(SD/sqrt(2),
    # sqrt(se_a^2 + se_b^2)/2); sigma = sqrt(se^2 + se_others^2).
    THREE_ONE_ODD = [(-8.584, 0.25), (-8.338, 0.285), (-7.284, 0.25)]

    def test_a_copy_a_decade_weaker_than_the_others_is_flagged(self):
        frame, nt, _ = _outlier_run({"A": self.THREE_ONE_ODD})
        # A3: delta 1.177; others' SE sqrt(0.0625 + 0.081225)/2 = 0.189555 beats 0.123.
        assert nt["A3"]["kdLog10SynonymDelta"] == pytest.approx(1.177, abs=1e-9)
        assert nt["A3"]["kdLog10SynonymDeltaSe"] == pytest.approx(0.313738, abs=1e-6)
        assert nt["A3"]["kdLog10SynonymZ"] == pytest.approx(3.7515432936671926, abs=1e-6)
        assert nt["A3"]["synonymOutlier"] is True
        # A1: delta -0.773; the others' spread 1.054/2 = 0.527 beats their fit errors.
        assert nt["A1"]["kdLog10SynonymDelta"] == pytest.approx(-0.773, abs=1e-9)
        assert nt["A1"]["kdLog10SynonymZ"] == pytest.approx(-1.3252378449411428, abs=1e-6)
        assert nt["A1"]["synonymOutlier"] is False
        assert nt["A2"]["kdLog10SynonymZ"] == pytest.approx(-0.5692258474565922, abs=1e-6)
        assert nt["A2"]["synonymOutlier"] is False
        assert _note(frame, "A") == "A3: weaker, +1.18 decades, z +3.75"

    def test_three_copies_that_agree_within_their_errors_are_not_flagged(self):
        frame, nt, summary = _outlier_run({"A": [(-6.177, 0.655), (-6.918, 0.27), (-6.551, 0.587)]})
        z = {key: row["kdLog10SynonymZ"] for key, row in nt.items()}
        assert z["A1"] == pytest.approx(0.7633461527395291, abs=1e-6)
        assert z["A2"] == pytest.approx(-1.073557733190834, abs=1e-6)
        assert z["A3"] == pytest.approx(-0.005042163883396568, abs=1e-6)
        assert max(abs(v) for v in z.values()) == pytest.approx(1.073557733190834, abs=1e-6)
        assert not any(row["synonymOutlier"] for row in nt.values())
        assert _note(frame, "A") is None
        assert summary["synonymOutlier"]["flaggedNtVariants"] == 0

    def test_two_copies_get_the_same_z_and_agree(self):
        frame, nt, _ = _outlier_run({"A": self.THREE_ONE_ODD[:2]})
        # The other copy's own SE: sigma = sqrt(0.25^2 + 0.285^2) = 0.379111 for both.
        assert nt["A1"]["kdLog10SynonymDeltaSe"] == pytest.approx(nt["A2"]["kdLog10SynonymDeltaSe"])
        assert nt["A1"]["kdLog10SynonymZ"] == pytest.approx(-0.6488868132152068, abs=1e-6)
        assert nt["A2"]["kdLog10SynonymZ"] == pytest.approx(0.6488868132152068, abs=1e-6)
        assert (nt["A1"]["synonymOutlier"], nt["A2"]["synonymOutlier"]) == (False, False)
        assert _note(frame, "A") is None

    def test_two_copies_that_disagree_are_named_as_a_pair(self):
        frame, nt, _ = _outlier_run({"A": [(-8.0, 0.1), (-7.0, 0.1)]})
        # |z| = 1 / sqrt(0.1^2 + 0.1^2) = 7.07 for both; neither is the odd one.
        assert (nt["A1"]["synonymOutlier"], nt["A2"]["synonymOutlier"]) == (True, True)
        assert nt["A2"]["kdLog10SynonymZ"] == pytest.approx(7.0710678118654755, abs=1e-6)
        assert _note(frame, "A") == "A1 vs A2: disagree, 1.00 decades apart, |z| 7.07"

    def test_without_a_se_the_dataset_sigma_scales_the_gap(self):
        groups = {
            "A": [(-8.0, None), (-8.2, None), (-8.1, None)],
            "B": [(-7.0, None), (-7.2, None)],
            "C": [(-9.0, None), (-9.1, None), (-9.2, None), (-7.9, None)],
        }
        frame, nt, summary = _outlier_run(groups)
        # Deviations from each mean, times sqrt(m / (m - 1)): A (0.1 -0.1 0) * sqrt(3/2);
        # B (0.1 -0.1) * sqrt(2); C (-0.2 -0.3 -0.4 0.9) * sqrt(4/3). Median -0.122474; absolute
        # deviations from it, median |C2 - center| = 0.223940: sigma 1.4826 * 0.223940 = 0.332007.
        sigma = 1.4826 * (0.3 * math.sqrt(4 / 3) - 0.1 * math.sqrt(3 / 2))
        assert summary["synonymOutlier"]["datasetSigma"] == pytest.approx(sigma, abs=1e-12)
        assert sigma == pytest.approx(0.33200703082780525, abs=1e-12)
        # C4: delta 1.2 over sigma * sqrt(1 + 1/3).
        assert nt["C4"]["kdLog10SynonymDelta"] == pytest.approx(1.2, abs=1e-9)
        assert nt["C4"]["kdLog10SynonymDeltaSe"] == pytest.approx(sigma * math.sqrt(4 / 3), abs=1e-12)
        assert nt["C4"]["kdLog10SynonymZ"] == pytest.approx(3.130146015131591, abs=1e-6)
        assert nt["C3"]["kdLog10SynonymZ"] == pytest.approx(-1.391176006725152, abs=1e-6)
        # B: delta 0.2 over sigma * sqrt(2).
        assert nt["B1"]["kdLog10SynonymZ"] == pytest.approx(0.4259589198599156, abs=1e-6)
        assert [key for key, row in nt.items() if row["synonymOutlier"]] == ["C4"]
        assert _note(frame, "C") == "C4: weaker, +1.20 decades, z +3.13"
        assert (_note(frame, "A"), _note(frame, "B")) == (None, None)
        assert summary["synonymOutlier"] == {
            "zThreshold": SYNONYM_OUTLIER_Z,
            "datasetSigma": pytest.approx(sigma, abs=1e-12),
            "testedAaVariants": 3,
            "testedNtVariants": 9,
            "flaggedNtVariants": 1,
            "flaggedAaVariants": 1,
        }

    def test_a_copy_at_a_wall_is_tested_at_the_wall_with_the_dataset_sigma(self):
        rows = [
            ("n1", 1e-8, "Good", False, 0.1),
            ("n2", 10.0**-8.1, "Good", False, 0.1),
            ("n3", 1e-13, "Good", True, None),
        ]
        links = pl.DataFrame({"ntKey": ["n1", "n2", "n3"], "aaKey": ["A", "A", "A"]})
        frame, per_nt, summary = combine_synonymous(
            _per_clonotype(rows), SynonymousInputs(links), log10_lo=LO, log10_hi=HI
        )
        nt = {r["ntKey"]: r for r in per_nt.iter_rows(named=True)}
        # Only n1 and n2 enter the sigma: deviations +-0.05 * sqrt(2), median 0, MAD 0.05 * sqrt(2).
        sigma = 1.4826 * 0.05 * math.sqrt(2)
        assert summary["synonymOutlier"]["datasetSigma"] == pytest.approx(sigma, abs=1e-12)
        # n3 enters at the lower wall, -10, has no SE, and is tested against n1 and n2.
        assert nt["n3"]["kdLog10SynonymDelta"] == pytest.approx(-1.95, abs=1e-9)
        assert nt["n3"]["kdLog10SynonymZ"] == pytest.approx(-1.95 / (sigma * math.sqrt(1.5)), abs=1e-9)
        # n1 is tested against n2 alone: the wall copy is never among the others.
        assert nt["n1"]["kdLog10SynonymDelta"] == pytest.approx(0.1, abs=1e-9)
        assert nt["n1"]["kdLog10SynonymZ"] == pytest.approx(0.1 / math.hypot(0.1, 0.1), abs=1e-9)
        assert _note(frame, "A") == "n3: tighter, at least 1.95 decades, below range"

    def test_a_copy_at_a_wall_does_not_drag_the_mean_the_others_are_tested_against(self):
        # No SEs: the dataset sigma throughout. Before the wall copy was dropped from the others,
        # A1 read against the mean of A2 and the wall, a decade away, and was flagged.
        groups = {"B": [(-7.0, None), (-7.2, None)], "C": [(-9.0, None), (-9.1, None), (-9.2, None)]}
        rows = [("A1", 1e-8, "Good", False, None), ("A2", 10.0**-8.1, "Good", False, None)]
        rows.append(("A3", 1e-13, "Good", True, None))
        links = [("A1", "A"), ("A2", "A"), ("A3", "A")]
        for aa, copies in groups.items():
            for i, (log10_kd, se) in enumerate(copies, start=1):
                rows.append((f"{aa}{i}", 10.0**log10_kd, "Good", False, se))
                links.append((f"{aa}{i}", aa))
        inputs = SynonymousInputs(pl.DataFrame(links, schema=["ntKey", "aaKey"], orient="row"))
        frame, per_nt, _ = combine_synonymous(_per_clonotype(rows), inputs, log10_lo=LO, log10_hi=HI)
        nt = {r["ntKey"]: r for r in per_nt.iter_rows(named=True)}
        assert nt["A1"]["kdLog10SynonymDelta"] == pytest.approx(0.1, abs=1e-9)
        assert (nt["A1"]["synonymOutlier"], nt["A2"]["synonymOutlier"], nt["A3"]["synonymOutlier"]) == (
            False,
            False,
            True,
        )
        assert _note(frame, "A") == "A3: tighter, at least 1.95 decades, below range"

    def test_a_measured_copy_and_a_copy_at_a_wall_test_the_wall_copy_alone(self):
        rows = [("n1", 1e-8, "Good", False, 0.1), ("n2", 1e-6 * 10, "Good", True, None)]
        rows += [("m1", 1e-8, "Good", False, 0.1), ("m2", 10.0**-8.2, "Good", False, 0.1)]
        links = pl.DataFrame({"ntKey": ["n1", "n2", "m1", "m2"], "aaKey": ["A", "A", "M", "M"]})
        frame, per_nt, _ = combine_synonymous(_per_clonotype(rows), SynonymousInputs(links), log10_lo=LO, log10_hi=HI)
        nt = {r["ntKey"]: r for r in per_nt.iter_rows(named=True)}
        # n1 has no other measured copy to read against.
        assert nt["n1"]["kdLog10SynonymDelta"] is None and nt["n1"]["synonymOutlier"] is None
        # n2 enters at the upper wall, -6: two decades above n1, at least.
        assert nt["n2"]["kdLog10SynonymDelta"] == pytest.approx(2.0, abs=1e-9)
        assert _note(frame, "A") == "n2: weaker, at least 2.00 decades, above range"

    def test_copies_at_a_wall_outnumbering_the_one_measured_copy_are_each_flagged_against_it(self):
        # Seen on Adams 2016: one copy measured, several at the upper wall. Each wall copy is tested
        # against the measured one alone and gets the same least gap; the measured copy has no other
        # measured copy and is not tested.
        rows = [("n1", 1e-8, "Good", False, 0.1)]
        rows += [(f"w{i}", 1e-5, "Above range", True, None) for i in (1, 2, 3)]
        rows += [("m1", 1e-8, "Good", False, 0.1), ("m2", 10.0**-8.2, "Good", False, 0.1)]
        links = pl.DataFrame({"ntKey": ["n1", "w1", "w2", "w3", "m1", "m2"], "aaKey": ["A", "A", "A", "A", "M", "M"]})
        frame, per_nt, summary = combine_synonymous(
            _per_clonotype(rows), SynonymousInputs(links), log10_lo=LO, log10_hi=HI
        )
        nt = {r["ntKey"]: r for r in per_nt.iter_rows(named=True)}
        # Only M's two copies enter the sigma: deviations +-0.1 * sqrt(2).
        sigma = 1.4826 * 0.1 * math.sqrt(2)
        assert summary["synonymOutlier"]["datasetSigma"] == pytest.approx(sigma, abs=1e-12)
        assert nt["n1"]["kdLog10SynonymDelta"] is None and nt["n1"]["synonymOutlier"] is None
        for w in ("w1", "w2", "w3"):
            assert nt[w]["kdLog10SynonymDelta"] == pytest.approx(2.0, abs=1e-9)
            assert nt[w]["kdLog10SynonymZ"] == pytest.approx(2.0 / (sigma * math.sqrt(2)), abs=1e-9)
            assert nt[w]["synonymOutlier"] is True
        assert _note(frame, "A") == "; ".join(f"w{i}: weaker, at least 2.00 decades, above range" for i in (1, 2, 3))
        # The variant is a mean, not a bound, and with its copies mostly at a wall it is Partial.
        a = frame.filter(pl.col("aaKey") == "A").row(0, named=True)
        assert (a["kdOutOfRange"], a["affinityClass"]) == (False, "Partial")
        assert summary["synonymOutlier"]["flaggedNtVariants"] == 3

    @pytest.mark.parametrize("k", [2, 3, 5])
    def test_the_dataset_sigma_recovers_the_spread_of_the_copies(self, k):
        # Each deviation is taken from a mean that includes the copy; unscaled, the MAD would read
        # sigma * sqrt((k - 1) / k): 0.71 at k = 2.
        rng = random.Random(7)
        groups = {
            f"A{g}": [(f"n{g}_{i}", Copy("measured", rng.gauss(-8.0, 0.3), None, "Good")) for i in range(k)]
            for g in range(3000)
        }
        assert dataset_sigma(groups) == pytest.approx(0.3, rel=0.05)

    def test_a_failed_copy_is_not_tested_and_does_not_count(self):
        rows = [("n1", 1e-8, "Good", False, 0.1), ("n2", 1e-7, "Good", False, 0.1), ("n3", None, "Failed", None, None)]
        links = pl.DataFrame({"ntKey": ["n1", "n2", "n3"], "aaKey": ["A", "A", "A"]})
        frame, per_nt, _ = combine_synonymous(_per_clonotype(rows), SynonymousInputs(links), log10_lo=LO, log10_hi=HI)
        nt = {r["ntKey"]: r for r in per_nt.iter_rows(named=True)}
        assert nt["n3"]["kdLog10SynonymDelta"] is None and nt["n3"]["synonymOutlier"] is None
        # k = 2: the pair.
        assert _note(frame, "A") == "n1 vs n2: disagree, 1.00 decades apart, |z| 7.07"

    def test_a_bound_a_single_copy_and_a_stop_are_not_tested(self):
        rows = [
            ("b1", 1e-13, "Good", True, None),
            ("b2", 1e-13, "Good", True, None),
            ("s1", 1e-8, "Good", False, 0.1),
            ("x1", 1e-8, "Good", False, 0.1),
            ("x2", 1e-7, "Good", False, 0.1),
        ]
        links = pl.DataFrame({"ntKey": ["b1", "b2", "s1", "x1", "x2"], "aaKey": ["B", "B", "S", "X", "X"]})
        productive = pl.DataFrame({"aaKey": ["B", "S", "X"], "productive": ["true", "true", "false"]})
        frame, per_nt, summary = combine_synonymous(
            _per_clonotype(rows), SynonymousInputs(links, aa_productive=productive), log10_lo=LO, log10_hi=HI
        )
        assert per_nt["kdLog10SynonymDelta"].null_count() == per_nt.height
        assert frame["synonymOutliers"].null_count() == frame.height
        assert summary["synonymOutlier"]["datasetSigma"] is None
        assert summary["synonymOutlier"]["testedAaVariants"] == 0

    def test_replicate_means_are_read_as_copies(self):
        # With replicates the copies are each nt variant's mean over its replicates.
        def merged(log10_kd, se, out_of_range=False):
            return Merged(log10_kd, se, None, Classification("Good", None), out_of_range, {})

        nt_values = {"n1": merged(-8.0, 0.1), "n2": merged(-8.1, 0.1), "n3": merged(-10.0, None, True)}
        aa_values = {"A": merged(-8.7, 0.6), "B": merged(-10.0, None, True)}
        mapping = {"n1": "A", "n2": "A", "n3": "A"}
        groups = outlier_groups(aa_values, nt_values, mapping, set())
        assert [nt for nt, _ in groups["A"]] == ["n1", "n2", "n3"]
        assert "B" not in groups
        # The out-of-range mean n3 is read as a wall copy: it stays out of the sigma.
        assert dataset_sigma(groups) == pytest.approx(1.4826 * 0.05 * math.sqrt(2), abs=1e-12)
        assert outlier_tests(groups["A"], dataset_sigma(groups))["n3"].bound is True

    def test_a_zero_sigma_gives_no_z(self):
        tests = outlier_tests(
            [("a", Copy("measured", -8.0, None, "Good")), ("b", Copy("measured", -8.0, None, "Good"))], 0.0
        )
        assert tests["a"] == OutlierTest(0.0, 0.0, None)
        assert tests["a"].flagged is None
        assert outlier_note(tests) is None
