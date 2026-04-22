"""Deterministic e2e corpus generator for titeseq-analysis.

Writes three reads parquets + a manifest describing the expected outcome for
every clonotype. Tests in tests/integration/test_corpus_e2e.py load the
artifacts and assert outputs against the manifest.

Run manually when the corpus needs to be regenerated:

    cd blocks/titeseq-analysis/software
    uv run python tests/fixtures/generate_corpus.py

Commit the regenerated parquets + manifest together with the code change
that required the update.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import polars as pl

CORPUS_DIR = Path(__file__).resolve().parent.parent / "data" / "corpus"
MASTER_SEED = 20260417

# Sub-µM dose grid (0.1 nM → 1 µM): realistic TiteSeq molar concentrations that
# stay below the attomolar-encoding ceiling enforced by R2 validation.
CONCENTRATIONS: list[float] = [0.0, 1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6]
BINS: list[int] = [1, 2, 3, 4]

# Shared Hill params for "nice" sigmoids — matches test_bin_mode_pipeline's noiseless-Hill defaults.
BASELINE = 1.5
AMPLITUDE = math.log(2.0)  # top ≈ baseline + 2.0 → 3.5 (near top of bin grid [1..4])
DEFAULT_SIGMA = 0.35  # Gaussian bin-assignment σ — matches pipeline regression tests.


def _rng(label: str) -> np.random.Generator:
    """Derive a per-clonotype rng so adding a new entry doesn't ripple-shift older ones."""
    seed = (MASTER_SEED + sum(ord(c) * (i + 1) for i, c in enumerate(label))) & 0xFFFFFFFF
    return np.random.default_rng(seed)


def _hill(x: float, baseline: float, amp: float, kd: float, n: float) -> float:
    top = baseline + math.exp(amp)
    if x <= 0:
        return baseline
    return baseline + (top - baseline) * (x**n) / (kd**n + x**n)


def _poisson_bin_reads(
    target_bin: float,
    per_conc: int,
    sigma: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Gaussian-weighted multinomial-ish allocation across BINS, Poisson-noised per bin."""
    centers = np.array(BINS, dtype=float)
    probs = np.exp(-0.5 * ((centers - target_bin) / sigma) ** 2)
    probs /= probs.sum()
    mean_reads = per_conc * probs
    return rng.poisson(mean_reads).astype(int)


def _bin_mode_rows(
    clonotype: str,
    targets: list[float],
    per_conc_list: list[int],
    sigma: float,
    rng: np.random.Generator,
    c0_sigma: float | None = None,
) -> list[dict]:
    rows: list[dict] = []
    for i, c in enumerate(CONCENTRATIONS):
        # Noisy clonotypes can still have clean c=0 so they don't skew the R6 global
        # baseline (arithmetic mean of c=0 signal across all clonotypes that survive
        # the R8 floor). Without this, wide-sigma clonotypes pull the baseline up
        # and depress the amplitude seen by every other clonotype's Hill fit.
        use_sigma = c0_sigma if c0_sigma is not None and c == 0.0 else sigma
        counts = _poisson_bin_reads(targets[i], per_conc_list[i], use_sigma, rng)
        for j, b in enumerate(BINS):
            rows.append(
                {
                    "clonotypeKey": clonotype,
                    "sampleId": f"{clonotype}_c{i}_b{b}",
                    "concentrationStr": str(c),
                    "concentration": float(c),
                    "bin": int(b),
                    "reads": int(counts[j]),
                }
            )
    return rows


def _sigmoid_targets(kd: float, n: float) -> list[float]:
    return [_hill(c, BASELINE, AMPLITUDE, kd, n) for c in CONCENTRATIONS]


def _expects_fit(entry: dict) -> bool:
    """True when the manifest entry anticipates a successful fit (Good or Partial outcome)."""
    if "expected_class" in entry:
        return entry["expected_class"] in ("Good", "Partial")
    if "expected_class_in" in entry:
        return any(c in ("Good", "Partial") for c in entry["expected_class_in"])
    return False


def _build_bin_mode_table() -> tuple[pl.DataFrame, dict]:
    """Assemble the multi-clonotype bin-mode reads table and its manifest entries."""
    entries: dict[str, dict] = {}
    all_rows: list[dict] = []

    n_conc = len(CONCENTRATIONS)

    # ---- G_LOW / G_MID / G_HIGH — clean Good cases across decades --------------
    # G_HIGH uses heavier reads so n converges tighter toward 1 (the decade-point
    # grid gives it fewer saturated points to anchor the plateau).
    for label, kd, per_conc, sig in [
        ("G_LOW", 1e-9, 500, DEFAULT_SIGMA),
        ("G_MID", 1e-8, 500, DEFAULT_SIGMA),
        ("G_HIGH", 1e-7, 3000, 0.25),  # denser + tighter σ; the kd=100 nM decade has fewer
        # saturated points on the grid so more signal + less noise is needed for a
        # clean Good classification.
    ]:
        rng = _rng(label)
        targets = _sigmoid_targets(kd, 1.0)
        all_rows += _bin_mode_rows(label, targets, [per_conc] * n_conc, sig, rng)
        entries[label] = {
            "expected_class": "Good",
            "expected_reason": None,
            "kd_range": [kd / 3.0, kd * 3.0],
            "hill_range": [0.5, 2.0],
            "kd_out_of_range": False,
            "hill_plot_position_is_sentinel": False,
            "notes": f"Clean sigmoid, kd ≈ {kd} — baseline Good across decades",
        }

    # ---- P_NOISY — wide σ bin assignment inflates residuals, R² should fall to Partial
    # c0_sigma=DEFAULT_SIGMA keeps the c=0 signal near BASELINE so this clonotype doesn't
    # skew the R6 global baseline (arithmetic mean across clonotypes) for every other fit.
    rng = _rng("P_NOISY")
    targets = _sigmoid_targets(1e-8, 1.0)
    all_rows += _bin_mode_rows(
        "P_NOISY", targets, [60] * n_conc, 0.85, rng, c0_sigma=DEFAULT_SIGMA
    )
    entries["P_NOISY"] = {
        "expected_class": "Partial",
        "expected_reason": None,
        "kd_range": [1e-9, 1e-7],
        "hill_range": [0.3, 10.0],
        "kd_out_of_range": False,
        "hill_plot_position_is_sentinel": False,
        "notes": "Wider σ bin assignment → weighted R² between r2_failed (0.5) and r2_good (0.8)",
    }

    # ---- F_LOW_R2 — clean sigmoid targets but heavy Gaussian bin-assignment noise.
    # Random scrambled targets risked triggering R9b (top-1 freq drop → hook); a
    # sigmoid + low per_conc + wide σ preserves the monotonic central trend while
    # blowing up per-concentration mean_bin noise, so the fit converges with R²
    # below r2_threshold_failed=0.5.
    rng = _rng("F_LOW_R2")
    targets = _sigmoid_targets(1e-8, 1.0)
    all_rows += _bin_mode_rows(
        "F_LOW_R2", targets, [200] * n_conc, 1.0, rng, c0_sigma=DEFAULT_SIGMA
    )
    entries["F_LOW_R2"] = {
        "expected_class": "Failed",
        "expected_reason_in": ["low_r2", "convergence_failure"],
        "kd_range": None,
        "hill_range": None,
        "kd_out_of_range": None,
        "hill_plot_position_is_sentinel": False,
        "notes": (
            "Sigmoid signal with noisy bin assignment → fit lands in Failed territory."
            " Accept either low_r2 or convergence_failure; both are valid expressions of"
            " 'signal too noisy for a useful fit'."
        ),
    }

    # ---- P_N_HIGH — steep Hill but high R² → Partial (reason None), n > n_max=2
    # (Truth-table row "r2 >= r2_good AND n out of range → Partial, reason None".)
    rng = _rng("P_N_HIGH")
    targets = _sigmoid_targets(1e-8, 4.0)
    all_rows += _bin_mode_rows("P_N_HIGH", targets, [500] * n_conc, DEFAULT_SIGMA, rng)
    entries["P_N_HIGH"] = {
        "expected_class": "Partial",
        "expected_reason": None,
        "kd_range": [3e-9, 3e-8],
        "hill_range": [2.0, 10.0],
        "kd_out_of_range": False,
        "hill_plot_position_is_sentinel": False,
        "notes": "Steep Hill (true n=4) + good data → high R², fit n > n_max=2 → Partial",
    }

    # ---- F_HOOK — R9b top-3 pattern, all top-3 reads well above hook_effect_min_reads
    # Top-1/2/3 concs are 1e-6/3e-7/1e-7 (indices 9/8/7 in CONCENTRATIONS).
    # Hit indices 7 & 8 with high signal (~3.5), drop index 9 to ~2.5 → top-2-top-1 = 1.0 > θ_bin=0.2.
    rng = _rng("F_HOOK")
    hook_targets = [
        1.5,  # c=0
        1.5,  # c=1e-10
        1.6,  # c=3e-10
        1.8,  # c=1e-9
        2.2,  # c=3e-9
        2.8,  # c=1e-8
        3.3,  # c=3e-8
        3.6,  # c=1e-7  (rank-3)
        3.6,  # c=3e-7  (rank-2)
        2.5,  # c=1e-6  (rank-1) ← dropped
    ]
    all_rows += _bin_mode_rows("F_HOOK", hook_targets, [500] * n_conc, DEFAULT_SIGMA, rng)
    entries["F_HOOK"] = {
        "expected_class": "Failed",
        "expected_reason": "non_monotonic_signal",
        "kd_range": None,
        "hill_range": None,
        "kd_out_of_range": None,
        "hill_plot_position_is_sentinel": True,
        "notes": "R9b: top-2 & top-3 elevated, top-1 drop > θ; every top-3 conc has reads ≥ hook_effect_min_reads",
    }

    # ---- F_INSUF_R — every concentration below read floor (all concs with 0 reads for this clonotype)
    # Construct as zero reads everywhere; clonotype is listed in the domain via a single non-zero cell at c=0
    rng = _rng("F_INSUF_R")
    # Put 1 read at each (c, b=1) so the clonotype exists in `all_clonotypes` but every conc fails floor.
    for i, c in enumerate(CONCENTRATIONS):
        for b in BINS:
            all_rows.append(
                {
                    "clonotypeKey": "F_INSUF_R",
                    "sampleId": f"F_INSUF_R_c{i}_b{b}",
                    "concentrationStr": str(c),
                    "concentration": float(c),
                    "bin": int(b),
                    "reads": 1 if b == 1 else 0,
                }
            )
    entries["F_INSUF_R"] = {
        "expected_class": "Failed",
        "expected_reason": "insufficient_reads",
        "kd_range": None,
        "hill_range": None,
        "kd_out_of_range": None,
        "hill_plot_position_is_sentinel": True,
        "notes": "All non-zero concs have clonotype_reads_at_conc < min_reads_per_concentration=3",
    }

    # ---- F_INSUF_P — only 3 non-zero concs pass the floor (< min_concentration_points=5)
    rng = _rng("F_INSUF_P")
    targets = _sigmoid_targets(1e-8, 1.0)
    # Heavily underfill most non-zero concs; only the top-3 concs survive the floor (3 < 5).
    # CONCENTRATIONS = [0, 1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6]
    per_conc = [200, 1, 1, 1, 1, 1, 1, 300, 300, 300]
    all_rows += _bin_mode_rows("F_INSUF_P", targets, per_conc, DEFAULT_SIGMA, rng)
    entries["F_INSUF_P"] = {
        "expected_class": "Failed",
        "expected_reason": "insufficient_points",
        "kd_range": None,
        "hill_range": None,
        "kd_out_of_range": None,
        "hill_plot_position_is_sentinel": True,
        "notes": "Only 3 non-zero concs pass read floor; below min_concentration_points=5",
    }

    # ---- F_FLAT — constant signal across all concentrations.
    # curve_fit converges (it can pick any kd + tiny amplitude with no data to contradict),
    # so kd/n are NOT null (not a sentinel case). R² is ~0 → Failed/low_r2. Included as
    # the "flat dose-response" dual of F_LOW_R2's "sigmoid + noisy".
    rng = _rng("F_FLAT")
    # Target c=0 at BASELINE so it doesn't pull the R6 global baseline; later
    # concentrations hold BASELINE too so there is no dose-response signal.
    flat_targets = [BASELINE] * n_conc
    all_rows += _bin_mode_rows("F_FLAT", flat_targets, [500] * n_conc, 0.08, rng)
    entries["F_FLAT"] = {
        "expected_class": "Failed",
        "expected_reason_in": ["low_r2", "convergence_failure"],
        "kd_range": None,
        "hill_range": None,
        "kd_out_of_range": None,
        "hill_plot_position_is_sentinel": False,
        "notes": (
            "Razor-flat signal anchored on BASELINE — no dose-response. Fit either"
            " fails to converge (amp_lo bound rejects flat signal) or converges at a"
            " rail kd with R² ≈ 0. Either Failed reason is acceptable."
        ),
    }

    # ---- K_LOW — Good-or-Partial fit, but kd below min non-zero concentration (0.1 nM)
    # kd=5e-11 (50 pM) is half the grid minimum — non-zero concs land in the plateau, so
    # the fit's n-coefficient is poorly constrained and often hits the N_HI rail. We accept
    # either Good or Partial; the point of this clonotype is R14b kdOutOfRange=true.
    rng = _rng("K_LOW")
    targets = _sigmoid_targets(5e-11, 1.0)
    all_rows += _bin_mode_rows("K_LOW", targets, [1500] * n_conc, DEFAULT_SIGMA, rng)
    entries["K_LOW"] = {
        "expected_class_in": ["Good", "Partial"],
        "expected_reason": None,
        "kd_range": [0.0, 1e-10],
        "hill_range": None,
        "kd_out_of_range": True,
        "hill_plot_position_is_sentinel": False,
        "notes": "Saturated Hill (true kd=5e-11) below min conc → kdOutOfRange=true; exact class depends on fit stability",
    }

    # K_HIGH (kd above max concentration) is not included in the corpus. The upper half
    # of R14b is covered by the integration/CLI hook test where null-kd Failed fits
    # land at the kdPlotPosition sentinel. K_LOW covers the lower bound of R14b here.

    # ---- C0_ONLY — reads only at c=0; no non-zero data points
    # Gaussian-distribute c=0 reads around BASELINE so this clonotype's c=0 signal
    # matches the global baseline and doesn't pull R6's mean upward.
    rng_c0 = _rng("C0_ONLY")
    c0_counts = _poisson_bin_reads(BASELINE, 500, DEFAULT_SIGMA, rng_c0)
    for j, b in enumerate(BINS):
        all_rows.append(
            {
                "clonotypeKey": "C0_ONLY",
                "sampleId": f"C0_ONLY_c0_b{b}",
                "concentrationStr": "0.0",
                "concentration": 0.0,
                "bin": int(b),
                "reads": int(c0_counts[j]),
            }
        )
    # Zero reads at every non-zero conc, so the clonotype exists in the all_clonotypes set.
    for i, c in enumerate(CONCENTRATIONS):
        if c == 0.0:
            continue
        for b in BINS:
            all_rows.append(
                {
                    "clonotypeKey": "C0_ONLY",
                    "sampleId": f"C0_ONLY_c{i}_b{b}",
                    "concentrationStr": str(c),
                    "concentration": float(c),
                    "bin": int(b),
                    "reads": 0,
                }
            )
    entries["C0_ONLY"] = {
        "expected_class": "Failed",
        "expected_reason": "insufficient_reads",
        "kd_range": None,
        "hill_range": None,
        "kd_out_of_range": None,
        "hill_plot_position_is_sentinel": True,
        "notes": "Only c=0 rows survive floor; zero non-zero points → insufficient_reads",
    }

    return pl.DataFrame(all_rows), entries


def _build_no_bin_mode_table() -> tuple[pl.DataFrame, dict]:
    """No-bin-mode corpus: 3 test clonotypes + a filler, all samples at non-zero
    concentrations only.

    Design: dropping c=0 means `compute_global_baseline` returns None and the pipeline
    uses the 4-parameter Hill fit (amplitude + baseline both free). Including c=0
    here would pollute the global arithmetic-mean baseline — one clonotype dominating
    a sample pulls the baseline off the per-clonotype library level and breaks every
    other fit. This mirrors the structure of the existing working no-bin integration
    test.
    """
    entries: dict[str, dict] = {}

    # Frequency profiles — indexed against NON_ZERO_CONCS below (not CONCENTRATIONS).
    non_zero_concs = [c for c in CONCENTRATIONS if c > 0]
    n_nz = len(non_zero_concs)

    # NB_GOOD: clean sigmoid ~0.01 → ~0.27 (kd near c=1e-8)
    good_freqs = [0.012, 0.018, 0.03, 0.05, 0.09, 0.15, 0.22, 0.25, 0.27]
    # NB_HOOK: rises then drops sharply at top conc → R9b hook flagged
    hook_freqs = [0.012, 0.018, 0.03, 0.05, 0.09, 0.15, 0.22, 0.25, 0.08]
    # NB_LOW_R2: scrambled — no consistent dose response
    lowr2_freqs = [0.05, 0.22, 0.03, 0.20, 0.04, 0.18, 0.05, 0.17, 0.06]
    assert all(len(f) == n_nz for f in (good_freqs, hook_freqs, lowr2_freqs))

    rng = np.random.default_rng(MASTER_SEED + 31)
    per_sample_depth = 20_000
    rows: list[dict] = []

    for i, c in enumerate(non_zero_concs):
        sample_id = f"no_bin_c{i + 1}"  # +1 keeps sample ids stable when c=0 re-added
        good_reads = int(rng.poisson(good_freqs[i] * per_sample_depth))
        hook_reads = int(rng.poisson(hook_freqs[i] * per_sample_depth))
        lowr2_reads = int(rng.poisson(lowr2_freqs[i] * per_sample_depth))
        # NB_FILLER absorbs the rest of the sample depth so the per-sample total
        # is stable and other clonotypes' frequencies reflect their intended profile.
        filler_reads = max(per_sample_depth - (good_reads + hook_reads + lowr2_reads), 0)

        for clonotype, rds in [
            ("NB_GOOD", good_reads),
            ("NB_HOOK", hook_reads),
            ("NB_LOW_R2", lowr2_reads),
            ("NB_FILLER", filler_reads),
        ]:
            rows.append(
                {
                    "clonotypeKey": clonotype,
                    "sampleId": sample_id,
                    "concentrationStr": str(c),
                    "concentration": float(c),
                    "reads": int(rds),
                }
            )

    entries["NB_GOOD"] = {
        "expected_class": "Good",
        "expected_reason": None,
        "kd_range": [1e-9, 2e-7],
        "hill_range": [0.3, 2.5],
        "kd_out_of_range": False,
        "hill_plot_position_is_sentinel": False,
        "notes": "R7b frequency signal — clean sigmoid rise to ~27% at max conc",
    }
    entries["NB_HOOK"] = {
        "expected_class": "Failed",
        "expected_reason": "non_monotonic_signal",
        "kd_range": None,
        "hill_range": None,
        "kd_out_of_range": None,
        "hill_plot_position_is_sentinel": True,
        "notes": "R9b no-bin threshold (0.02): top-1 freq < top-2,top-3 → hook flagged",
    }
    entries["NB_LOW_R2"] = {
        "expected_class": "Failed",
        "expected_reason_in": ["low_r2", "convergence_failure"],
        "kd_range": None,
        "hill_range": None,
        "kd_out_of_range": None,
        "hill_plot_position_is_sentinel": False,
        "notes": (
            "Scrambled frequencies → Hill fit either fails to converge or converges to"
            " low R². Either Failed reason is acceptable; the branch under test is"
            " 'scrambled no-bin signal does not classify as Good or Partial'."
            " Sentinel positions only apply on convergence_failure, so not asserted."
        ),
    }
    # NB_FILLER is not in the manifest expectations — it exists only to stabilize
    # per-sample depth so test clonotype frequencies reflect their intended profile.

    return pl.DataFrame(rows), entries


def _build_antigen_table() -> tuple[pl.DataFrame, dict]:
    """R4 antigen filter corpus: two target clonotypes, two distractors to be filtered out."""
    entries: dict[str, dict] = {}
    rows: list[dict] = []

    def append(clonotype: str, antigen: str, kd: float) -> None:
        rng = _rng(f"ANT_{clonotype}")
        targets = _sigmoid_targets(kd, 1.0)
        batch = _bin_mode_rows(clonotype, targets, [500] * len(CONCENTRATIONS), 0.25, rng)
        for r in batch:
            r["antigen"] = antigen
        rows.extend(batch)

    append("ANT_T1", "target", 5e-9)
    append("ANT_T2", "target", 5e-8)
    append("ANT_D1", "other", 5e-9)
    append("ANT_D2", "other", 5e-8)

    entries["ANT_T1"] = {
        "expected_class": "Good",
        "expected_reason": None,
        "kd_range": [1.5e-9, 1.5e-8],
        "hill_range": [0.5, 2.0],
        "kd_out_of_range": False,
        "hill_plot_position_is_sentinel": False,
        "notes": "antigen=target, kd≈5 nM — survives R4 filter",
    }
    entries["ANT_T2"] = {
        "expected_class": "Good",
        "expected_reason": None,
        "kd_range": [1.5e-8, 1.5e-7],
        "hill_range": [0.5, 2.0],
        "kd_out_of_range": False,
        "hill_plot_position_is_sentinel": False,
        "notes": "antigen=target, kd≈50 nM — survives R4 filter",
    }
    # Distractors have no manifest entry beyond "must not appear in output".
    return pl.DataFrame(rows), entries


def main() -> None:
    CORPUS_DIR.mkdir(parents=True, exist_ok=True)

    bin_df, bin_entries = _build_bin_mode_table()
    no_bin_df, no_bin_entries = _build_no_bin_mode_table()
    antigen_df, antigen_entries = _build_antigen_table()

    bin_df.write_parquet(CORPUS_DIR / "reads_bin_mode.parquet")
    no_bin_df.write_parquet(CORPUS_DIR / "reads_no_bin_mode.parquet")
    antigen_df.write_parquet(CORPUS_DIR / "reads_antigen.parquet")

    manifest = {
        "corpus_version": "1.0.0",
        "master_seed": MASTER_SEED,
        "concentrations": CONCENTRATIONS,
        "bins": BINS,
        "max_non_zero_concentration": max(c for c in CONCENTRATIONS if c > 0),
        "bin_mode": {
            "clonotypes": bin_entries,
            "fitted_mean_bin_clonotypes": sorted(
                name for name, entry in bin_entries.items() if _expects_fit(entry)
            ),
        },
        "no_bin_mode": {
            "clonotypes": no_bin_entries,
            "fitted_mean_bin_clonotypes": sorted(
                name for name, entry in no_bin_entries.items() if _expects_fit(entry)
            ),
        },
        "antigen": {
            "target_antigen": "target",
            "clonotypes": antigen_entries,
            "distractor_clonotypes": ["ANT_D1", "ANT_D2"],
        },
    }
    (CORPUS_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    print(f"Wrote corpus to {CORPUS_DIR}")
    print(f"  bin_mode:    {bin_df.height} rows, {bin_df['clonotypeKey'].n_unique()} clonotypes")
    print(f"  no_bin_mode: {no_bin_df.height} rows, {no_bin_df['clonotypeKey'].n_unique()} clonotypes")
    print(f"  antigen:     {antigen_df.height} rows, {antigen_df['clonotypeKey'].n_unique()} clonotypes")


if __name__ == "__main__":
    main()
