"""The parent delta of a run with sample replicates: taken inside each replicate, then combined."""

from __future__ import annotations

import math

import polars as pl
import pytest

from replicate_combine import ReplicateValue, combine_replicates
from synonymous import COL_AA, COL_MUTATIONS, COL_PARENT, parent_deltas

LO, HI = -9.5, -5.0
REPLICATES = ("r1", "r2", "r3")

LINKS = pl.DataFrame({COL_AA: ["P", "V"], COL_PARENT: ["p0", "p0"]})
MUTATIONS = pl.DataFrame({COL_AA: ["P", "V"], COL_MUTATIONS: [0, 1]})


def measured(replicate: str, log10_kd: float, se: float | None = 0.1) -> ReplicateValue:
    return ReplicateValue("measured", log10_kd, se, "Good", replicate=replicate)


def failed(replicate: str) -> ReplicateValue:
    return ReplicateValue("failed", None, None, "Failed", replicate=replicate)


def above(replicate: str) -> ReplicateValue:
    return ReplicateValue("above", HI, None, "Above range", replicate=replicate)


def deltas(parent: list[ReplicateValue], variant: list[ReplicateValue]):
    """(delta, SE) of the parent and of the variant, and the counts."""
    replicates = {"P": parent, "V": variant}
    values = {key: combine_replicates(members) for key, members in replicates.items()}
    columns, counts = parent_deltas(values, LINKS, MUTATIONS, COL_AA, replicates=replicates)
    d, se = columns["kdLog10ParentDelta"], columns["kdLog10ParentDeltaSe"]
    return (d[0], se[0]), (d[1], se[1]), counts


# The parent's value moves by replicate; the variant sits one decade above it in each.
PARENT = [measured(r, v) for r, v in zip(REPLICATES, (-8.0, -7.5, -8.4))]
VARIANT = [measured(r, v + 1.0) for r, v in zip(REPLICATES, (-8.0, -7.5, -8.4))]


def test_a_shift_shared_by_a_replicates_variant_and_parent_cancels():
    _, (delta, se), _ = deltas(PARENT, VARIANT)
    assert delta == pytest.approx(1.0)
    # No spread is left, so the SE is the fit-error term: RMS of hypot(0.1, 0.1), over sqrt(3).
    assert se == pytest.approx(math.hypot(0.1, 0.1) / math.sqrt(3))


def test_a_replicate_where_the_variant_failed_is_left_out_of_the_delta():
    _, (delta, se), _ = deltas(PARENT, [VARIANT[0], failed("r2"), VARIANT[2]])
    assert delta == pytest.approx(1.0)
    # Two replicates are left, and the divisor counts those two.
    assert se == pytest.approx(math.hypot(0.1, 0.1) / math.sqrt(2))


@pytest.mark.parametrize(
    ("parent", "variant"),
    [
        ([PARENT[0], above("r2"), PARENT[2]], VARIANT),
        (PARENT, [VARIANT[0], VARIANT[1], above("r3")]),
    ],
    ids=["parent at a range limit", "variant at a range limit"],
)
def test_a_replicate_at_a_range_limit_sends_the_delta_to_the_combined_values(parent, variant):
    # That replicate is in the combined value at the limit. A delta from the other replicates alone
    # could contradict the two values it stands beside.
    _, (delta, se), counts = deltas(parent, variant)
    combined = combine_replicates(variant).log10_kd - combine_replicates(parent).log10_kd
    assert delta == pytest.approx(combined)
    assert se is None
    assert (counts["replicateDeltas"], counts["combinedValueDeltas"]) == (0, 1)


def test_a_variant_tighter_than_its_parent_on_the_combined_values_has_a_negative_delta():
    # One replicate measures it a decade weaker; the other two put it under the lower limit.
    parent = [measured(r, -8.0) for r in REPLICATES]
    below = [ReplicateValue("below", LO, None, "Below range", replicate=r) for r in ("r2", "r3")]
    _, (delta, _), _ = deltas(parent, [measured("r1", -7.0), *below])
    assert delta == pytest.approx((-7.0 + 2 * LO) / 3 - -8.0)
    assert delta < 0


def test_fits_without_an_se_give_the_spread_of_the_deltas():
    parent = [measured(r, v, se=None) for r, v in zip(REPLICATES, (-8.0, -7.5, -8.4))]
    variant = [measured(r, v, se=None) for r, v in zip(REPLICATES, (-7.0, -6.0, -7.9))]
    _, (delta, se), _ = deltas(parent, variant)
    assert delta == pytest.approx(1.0)
    assert se == pytest.approx(0.5 / math.sqrt(3))


def test_replicates_missing_a_fit_se_count_at_the_typical_error_of_the_others():
    # Deltas all 1.0, so the fit-error term rules. r1 and r2 report hypot(0.1, 0.1) and
    # hypot(0.3, 0.1); r3 reports none. RMS over the two that report, over sqrt(3).
    parent = [measured("r1", -8.0), measured("r2", -7.5), measured("r3", -8.4, se=None)]
    variant = [measured("r1", -7.0), measured("r2", -6.5, se=0.3), measured("r3", -7.4)]
    _, (delta, se), _ = deltas(parent, variant)
    assert delta == pytest.approx(1.0)
    assert se == pytest.approx(math.sqrt((0.02 + 0.10) / 2) / math.sqrt(3))


def test_one_shared_replicate_without_a_fit_se_has_no_se():
    parent = [measured("r1", -8.0, se=None), failed("r2"), failed("r3")]
    variant = [measured(r, v, se=None) for r, v in zip(REPLICATES, (-7.0, -6.5, -7.4))]
    _, (delta, se), _ = deltas(parent, variant)
    assert delta == pytest.approx(1.0)
    assert se is None


def test_the_parents_own_value_stands_beside_the_delta():
    replicates = {"P": PARENT, "V": VARIANT}
    values = {key: combine_replicates(members) for key, members in replicates.items()}
    columns, _ = parent_deltas(values, LINKS, MUTATIONS, COL_AA, replicates=replicates)
    assert columns["kdParent"][1] == pytest.approx(10.0 ** values["P"].log10_kd)
    assert columns["kdLog10ParentSe"][1] == pytest.approx(values["P"].log10_kd_se)


def test_the_spread_of_the_deltas_sets_the_se_when_it_is_the_larger_term():
    variant = [measured(r, v) for r, v in zip(REPLICATES, (-7.0, -6.0, -7.9))]
    _, (delta, se), _ = deltas(PARENT, variant)
    per_replicate = [1.0, 1.5, 0.5]
    assert delta == pytest.approx(1.0)
    sd = math.sqrt(sum((d - 1.0) ** 2 for d in per_replicate) / 2)
    assert se == pytest.approx(sd / math.sqrt(3))


def test_one_shared_replicate_reports_its_own_delta_and_fit_error():
    _, (delta, se), _ = deltas([PARENT[0], failed("r2"), failed("r3")], VARIANT)
    assert delta == pytest.approx(1.0)
    assert se == pytest.approx(math.hypot(0.1, 0.1))


def test_with_no_shared_measured_replicate_the_combined_values_give_a_delta_without_an_se():
    parent = [PARENT[0], failed("r2"), failed("r3")]
    variant = [failed("r1"), VARIANT[1], VARIANT[2]]
    _, (delta, se), counts = deltas(parent, variant)
    assert delta == pytest.approx((-6.5 + -7.4) / 2 - -8.0)
    assert se is None
    assert counts["combinedValueDeltas"] == 1


def test_the_parent_itself_is_zero():
    (delta, se), _, _ = deltas(PARENT, VARIANT)
    assert (delta, se) == (0.0, 0.0)


def test_a_variant_with_no_combined_value_has_no_delta():
    _, (delta, se), _ = deltas(PARENT, [failed(r) for r in REPLICATES])
    assert (delta, se) == (None, None)


def test_the_counts_say_how_many_deltas_came_from_the_replicates():
    _, _, counts = deltas(PARENT, VARIANT)
    assert counts["replicateDeltas"] == 1
    assert counts["combinedValueDeltas"] == 0


def test_without_replicates_the_delta_is_the_difference_of_the_values():
    values = {"P": combine_replicates(PARENT), "V": combine_replicates(VARIANT)}
    columns, counts = parent_deltas(values, LINKS, MUTATIONS, COL_AA)
    assert columns["kdLog10ParentDelta"][1] == pytest.approx(1.0)
    assert columns["kdLog10ParentDeltaSe"][1] == pytest.approx(
        math.hypot(values["P"].log10_kd_se, values["V"].log10_kd_se)
    )
    assert "replicateDeltas" not in counts
