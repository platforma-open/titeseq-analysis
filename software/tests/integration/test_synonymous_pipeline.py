"""Synonymous nucleotide copies through the whole pipeline, on a synthetic library with known K_D.

Each amino-acid variant is encoded by n nucleotide copies whose reads differ in depth and spread, so
their fits differ. Expectations are bands around the true K_D, and exact counts.
"""

from __future__ import annotations

import json
import math

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from constants import FitParams
from main import main
from pipeline import run
from synonymous import SynonymousInputs

CONCS = [1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6]
BINS = [1, 2, 3, 4]
GATE_WIDTH = 0.65
LOG_LO = math.log10(min(CONCS))
LOG_HI = math.log10(max(CONCS))

LF = FitParams(equation="log_fluorescence")
HILL = FitParams(equation="hill_gate_index")


def _reads(nt: str, true_kd: float, *, per_conc: int = 800, sd_bins: float = 0.35) -> list[dict]:
    """Reads drawn from log10(A*c/(c+K) + B), spread over the gates as a Gaussian in gate index."""
    span = (len(BINS) - 1) * GATE_WIDTH
    saturated, baseline = 10.0 ** (0.85 * span), 2.0
    c = np.asarray(CONCS, dtype=float)
    targets = np.log10((saturated - baseline) * c / (c + true_kd) + baseline) / GATE_WIDTH + BINS[0]
    rows = []
    for i, conc in enumerate(CONCS):
        weights = np.exp(-0.5 * ((np.array(BINS, dtype=float) - targets[i]) / sd_bins) ** 2)
        counts = np.round(weights / weights.sum() * per_conc).astype(int)
        for j, b in enumerate(BINS):
            rows.append(
                {
                    "clonotypeKey": nt,
                    "sampleId": f"s_c{i}_b{b}",
                    "concentrationStr": str(conc),
                    "concentration": float(conc),
                    "bin": int(b),
                    "reads": int(counts[j]),
                }
            )
    return rows


# aa variant -> [(nt copy, true K_D, depth, spread)]
LIBRARY = {
    "WT": [("wt_a", 1e-8, 800, 0.35), ("wt_b", 1e-8, 600, 0.40), ("wt_c", 1e-8, 1000, 0.30)],
    "TIGHT": [("t_a", 1e-14, 800, 0.35), ("t_b", 1e-14, 600, 0.40)],
    "DEAD": [("d_a", 1e3, 800, 0.35), ("d_b", 1e3, 600, 0.40)],
    "MIXED": [("x_a", 1e-8, 800, 0.35), ("x_b", 1e-14, 800, 0.35)],
    "SINGLE": [("s_a", 1e-9, 800, 0.35)],
    "LOST": [("l_a", 1e-8, 2, 0.35)],
}
# Fitted, but not the copy of exactly one aa variant.
STRAYS = [("amb", 1e-8, 800, 0.35), ("orphan", 1e-8, 800, 0.35)]


def _library_reads() -> pl.DataFrame:
    copies = [c for cs in LIBRARY.values() for c in cs] + STRAYS
    return pl.DataFrame([row for nt, kd, depth, sd in copies for row in _reads(nt, kd, per_conc=depth, sd_bins=sd)])


def _links() -> pl.DataFrame:
    pairs = [(nt, aa) for aa, cs in LIBRARY.items() for nt, *_ in cs] + [("amb", "WT"), ("amb", "SINGLE")]
    return pl.DataFrame({"ntKey": [p[0] for p in pairs], "aaKey": [p[1] for p in pairs]})


PARENT_INPUTS = {
    "aa_mutations": pl.DataFrame({"aaKey": list(LIBRARY), "mutationCount": ["0", "1", "1", "1", "1", "1"]}),
    "aa_parent_links": pl.DataFrame({"aaKey": list(LIBRARY), "parentId": ["P1"] * len(LIBRARY)}),
}


@pytest.fixture(scope="module")
def reads() -> pl.DataFrame:
    return _library_reads()


@pytest.fixture(scope="module")
def lf(reads):
    return run(reads, params=LF, synonymous=SynonymousInputs(_links(), **PARENT_INPUTS))


def _aa(outputs, aa: str) -> dict:
    return outputs["per_variant"].filter(pl.col("aaKey") == aa).row(0, named=True)


class TestEachCopyIsFittedOnItsOwn:
    @pytest.mark.parametrize("params", [LF, HILL], ids=["log_fluorescence", "hill_gate_index"])
    def test_the_nucleotide_outputs_are_those_of_a_run_without_the_linker(self, reads, params):
        plain = run(reads, params=params)
        combined = run(reads, params=params, synonymous=SynonymousInputs(_links(), **PARENT_INPUTS))
        for name in ("per_clonotype", "mean_bin", "fitted_mean_bin", "concentration_value"):
            assert_frame_equal(combined[name], plain[name])
        assert set(combined) - set(plain) == {"per_variant", "synonymous_summary"}

    def test_the_aa_value_is_the_mean_of_its_copies_fits_not_a_fit_of_pooled_reads(self, lf):
        per_nt = lf["per_clonotype"].filter(pl.col("clonotypeKey").is_in(["wt_a", "wt_b", "wt_c"]))
        mean = sum(math.log10(kd) for kd in per_nt["kd"].to_list()) / 3
        assert math.log10(_aa(lf, "WT")["kd"]) == pytest.approx(mean, abs=1e-9)


class TestLogFluorescence:
    def test_three_measured_copies_recover_the_kd_with_a_scatter_se(self, lf):
        wt = _aa(lf, "WT")
        assert math.log10(wt["kd"]) == pytest.approx(-8.0, abs=0.3)
        assert (wt["ntVariants"], wt["ntMeasured"], wt["ntFailed"]) == (3, 3, 0)
        assert 0.0 <= wt["kdLog10Se"] <= 0.5
        assert wt["affinityClass"] == "Good"
        assert wt["kdOutOfRange"] is False

    def test_the_se_is_the_larger_of_the_spread_and_the_copies_fit_errors(self, lf):
        copies = lf["per_clonotype"].filter(pl.col("clonotypeKey").is_in(["wt_a", "wt_b", "wt_c"]))
        values = [math.log10(kd) for kd in copies["kd"].to_list()]
        spread = float(np.std(values, ddof=1)) / math.sqrt(3)
        fit_term = math.sqrt(sum(se**2 for se in copies["kdLog10Se"].to_list())) / 3
        assert _aa(lf, "WT")["kdLog10Se"] == pytest.approx(max(spread, fit_term), rel=1e-9)

    def test_copies_saturated_across_the_titration_give_a_lower_bound(self, lf):
        tight = _aa(lf, "TIGHT")
        assert math.log10(tight["kd"]) == pytest.approx(LOG_LO, abs=1e-9)
        assert (tight["ntCensoredBelow"], tight["ntCensoredAbove"]) == (2, 0)
        assert tight["kdLog10Se"] is None
        assert tight["kdOutOfRange"] is True
        assert tight["affinityClass"] == "Partial"

    def test_flat_dim_copies_give_an_upper_bound(self, lf):
        dead = _aa(lf, "DEAD")
        assert math.log10(dead["kd"]) == pytest.approx(LOG_HI, abs=1e-9)
        assert (dead["ntCensoredBelow"], dead["ntCensoredAbove"]) == (0, 2)
        assert dead["affinityClass"] == "Partial"

    def test_a_censored_copy_enters_at_its_wall_beside_a_measured_one(self, lf):
        mixed = _aa(lf, "MIXED")
        assert (mixed["ntMeasured"], mixed["ntCensoredBelow"]) == (1, 1)
        measured = lf["per_clonotype"].filter(pl.col("clonotypeKey") == "x_a")["kd"].item()
        assert math.log10(mixed["kd"]) == pytest.approx((math.log10(measured) + LOG_LO) / 2, abs=1e-9)
        assert mixed["kdOutOfRange"] is False

    def test_one_copy_carries_its_own_fit_se(self, lf):
        single = _aa(lf, "SINGLE")
        own = lf["per_clonotype"].filter(pl.col("clonotypeKey") == "s_a").row(0, named=True)
        assert single["ntVariants"] == 1
        assert single["kdLog10Se"] == own["kdLog10Se"]
        assert single["kd"] == pytest.approx(own["kd"], rel=1e-12)
        assert math.log10(single["kd"]) == pytest.approx(-9.0, abs=0.3)

    def test_a_variant_whose_only_copy_failed_is_failed_with_no_value(self, lf):
        lost = _aa(lf, "LOST")
        assert (lost["ntVariants"], lost["ntFailed"]) == (1, 1)
        assert lost["kd"] is None
        assert lost["affinityClass"] == "Failed"

    def test_unmapped_and_ambiguous_copies_are_counted_not_averaged(self, lf):
        summary = lf["synonymous_summary"]
        assert summary["ntVariants"] == 13
        assert summary["aaVariants"] == len(LIBRARY)
        assert (summary["unmappedNtVariants"], summary["ambiguousNtVariants"]) == (1, 1)
        assert _aa(lf, "WT")["ntVariants"] == 3


class TestParentDelta:
    def test_the_parent_row_is_zero_with_zero_error(self, lf):
        wt = _aa(lf, "WT")
        assert (wt["kdLog10ParentDelta"], wt["kdLog10ParentDeltaSe"]) == (0.0, 0.0)

    def test_a_variant_reads_against_its_parent(self, lf):
        single, wt = _aa(lf, "SINGLE"), _aa(lf, "WT")
        assert single["kdLog10ParentDelta"] == pytest.approx(math.log10(single["kd"]) - math.log10(wt["kd"]), abs=1e-9)
        assert single["kdLog10ParentDeltaSe"] == pytest.approx(
            math.hypot(single["kdLog10Se"], wt["kdLog10Se"]), abs=1e-12
        )
        assert single["kdLog10ParentDelta"] == pytest.approx(-1.0, abs=0.4)

    def test_a_failed_variant_has_no_delta(self, lf):
        assert _aa(lf, "LOST")["kdLog10ParentDelta"] is None

    def test_without_parent_columns_there_are_no_delta_columns(self, reads):
        out = run(reads, params=LF, synonymous=SynonymousInputs(_links()))
        assert "kdLog10ParentDelta" not in out["per_variant"].columns
        assert out["synonymous_summary"]["parentMode"] is False


class TestReadWeighting:
    @pytest.fixture(scope="class")
    def weighted(self, reads):
        return run(reads, params=LF, synonymous=SynonymousInputs(_links(), weighting="reads", **PARENT_INPUTS))

    def test_each_copy_weighs_its_total_reads_over_the_titration(self, reads, weighted):
        totals = {nt: reads.filter(pl.col("clonotypeKey") == nt)["reads"].sum() for nt in ("wt_a", "wt_b", "wt_c")}
        fits = weighted["per_clonotype"].filter(pl.col("clonotypeKey").is_in(list(totals)))
        x = {r["clonotypeKey"]: math.log10(r["kd"]) for r in fits.iter_rows(named=True)}
        mean = sum(totals[nt] * x[nt] for nt in totals) / sum(totals.values())
        assert math.log10(_aa(weighted, "WT")["kd"]) == pytest.approx(mean, abs=1e-9)
        assert weighted["synonymous_summary"]["weighting"] == "reads"

    def test_the_weighted_mean_differs_from_the_plain_mean(self, lf, weighted):
        assert _aa(weighted, "WT")["kd"] != pytest.approx(_aa(lf, "WT")["kd"], rel=1e-12)

    def test_the_nucleotide_outputs_do_not_move(self, lf, weighted):
        assert_frame_equal(weighted["per_clonotype"], lf["per_clonotype"])


class TestHillGateIndex:
    def test_one_copy_has_no_se_and_keeps_its_own_class(self, reads):
        out = run(reads, params=HILL, synonymous=SynonymousInputs(_links()))
        single = _aa(out, "SINGLE")
        own = out["per_clonotype"].filter(pl.col("clonotypeKey") == "s_a").row(0, named=True)
        assert single["kdLog10Se"] is None
        assert single["affinityClass"] == own["affinityClass"]

    def test_several_copies_get_a_scatter_se(self, reads):
        out = run(reads, params=HILL, synonymous=SynonymousInputs(_links()))
        wt = _aa(out, "WT")
        copies = out["per_clonotype"].filter(pl.col("clonotypeKey").is_in(["wt_a", "wt_b", "wt_c"]))
        assert copies["kdOutOfRange"].to_list() == [False, False, False]
        assert wt["ntMeasured"] == 3
        assert wt["kdLog10Se"] is not None
        mean = sum(math.log10(kd) for kd in copies["kd"].to_list()) / 3
        assert math.log10(wt["kd"]) == pytest.approx(mean, abs=1e-9)


def test_the_cli_writes_the_amino_acid_frame_and_the_summary(tmp_path, reads):
    reads.write_parquet(tmp_path / "reads.parquet")
    _links().write_csv(tmp_path / "nt_to_aa.tsv", separator="\t")
    PARENT_INPUTS["aa_mutations"].write_csv(tmp_path / "aa_mutations.tsv", separator="\t")
    PARENT_INPUTS["aa_parent_links"].write_csv(tmp_path / "aa_parent_links.tsv", separator="\t")
    (tmp_path / "params.json").write_text(json.dumps({"equation": "log_fluorescence", "plateau_floor": 2 / 3}))
    out = {name: str(tmp_path / f"{name}.tsv") for name in ("pc", "mb", "fmb", "cv", "pv")}
    rc = main(
        [
            "--reads", str(tmp_path / "reads.parquet"),
            "--params", str(tmp_path / "params.json"),
            "--out-per-clonotype", out["pc"],
            "--out-mean-bin", out["mb"],
            "--out-fitted-mean-bin", out["fmb"],
            "--out-concentration-value", out["cv"],
            "--nt-to-aa", str(tmp_path / "nt_to_aa.tsv"),
            "--aa-mutations", str(tmp_path / "aa_mutations.tsv"),
            "--aa-parent-links", str(tmp_path / "aa_parent_links.tsv"),
            "--out-per-variant", out["pv"],
            "--out-synonymous-summary", str(tmp_path / "summary.json"),
        ]
    )  # fmt: skip
    assert rc == 0
    per_variant = pl.read_csv(out["pv"], separator="\t")
    assert per_variant["aaKey"].to_list() == sorted(LIBRARY)
    assert "kdLog10ParentDeltaSe" in per_variant.columns
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["parentMode"] is True
    assert summary["unidentifiableParents"] == []


def test_the_cli_passes_read_weighting_on(tmp_path, reads):
    reads.write_parquet(tmp_path / "reads.parquet")
    _links().write_csv(tmp_path / "nt_to_aa.tsv", separator="\t")
    (tmp_path / "params.json").write_text(json.dumps({"equation": "log_fluorescence", "plateau_floor": 2 / 3}))
    out = {name: str(tmp_path / f"{name}.tsv") for name in ("pc", "mb", "fmb", "cv", "pv")}
    rc = main(
        [
            "--reads", str(tmp_path / "reads.parquet"),
            "--params", str(tmp_path / "params.json"),
            "--out-per-clonotype", out["pc"],
            "--out-mean-bin", out["mb"],
            "--out-fitted-mean-bin", out["fmb"],
            "--out-concentration-value", out["cv"],
            "--nt-to-aa", str(tmp_path / "nt_to_aa.tsv"),
            "--out-per-variant", out["pv"],
            "--out-synonymous-summary", str(tmp_path / "summary.json"),
            "--synonymous-weighting", "reads",
        ]
    )  # fmt: skip
    assert rc == 0
    assert json.loads((tmp_path / "summary.json").read_text())["weighting"] == "reads"


def test_the_cli_refuses_weighting_without_a_linker():
    with pytest.raises(SystemExit):
        main(
            [
                "--reads", "r.parquet",
                "--out-per-clonotype", "a", "--out-mean-bin", "b",
                "--out-fitted-mean-bin", "c", "--out-concentration-value", "d",
                "--synonymous-weighting", "reads",
            ]
        )  # fmt: skip


def test_the_cli_refuses_a_linker_without_its_outputs(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "--reads", "r.parquet",
                "--out-per-clonotype", "a", "--out-mean-bin", "b",
                "--out-fitted-mean-bin", "c", "--out-concentration-value", "d",
                "--nt-to-aa", "links.tsv",
            ]
        )  # fmt: skip
