# Deterministic E2E Corpus

Committed reads tables + per-clonotype manifest used by `tests/integration/test_corpus_e2e.py` to exercise the full titeseq pipeline (`io_layer → normalization → pre_fit → hill_fit → classify → output_build`) against known-truth inputs.

## Files

| File | Purpose |
|---|---|
| `reads_bin_mode.parquet` | 12 clonotypes covering every classification branch in bin mode |
| `reads_no_bin_mode.parquet` | 4 clonotypes covering R7b frequency signal + NB-mode hook detection |
| `reads_antigen.parquet` | 4 clonotypes exercising the R4 antigen filter (2 targets + 2 distractors) |
| `manifest.json` | Per-clonotype expected outcome (class, reason, kd/n tolerance bands, flags) |

## Bin-mode clonotypes

| Key | Intent |
|---|---|
| `G_LOW` / `G_MID` / `G_HIGH` | Clean sigmoid Good cases at kd ≈ 1 / 10 / 100 nM |
| `P_NOISY` | Wide-σ bin assignment → weighted R² in the Partial band |
| `P_N_HIGH` | Steep true Hill (n = 4) + good data. Default run: n pinned at 1 → Good. With `pin_hill_exponent=False`: n > `n_max` → Partial |
| `K_LOW` | Saturated Hill below min non-zero conc → clamped to that conc, `CensoredBelow`, `kdCensor = "below"` |
| `F_LOW_R2` | Sigmoid + heavy bin-assignment noise → R² below `r2_threshold_failed`, still reported. K_D rails to the bottom → `CensoredBelow` |
| `F_FLAT` | Razor-flat signal → the fit lands on the amplitude floor → Failed/insufficient_dynamic_range |
| `F_HOOK` | R9b pattern: top-1 drops, top-2/3 elevated. Default run: detector off → Partial. With `hook_detector_enabled=True`: Failed/non_monotonic_signal |
| `F_INSUF_R` | Every non-zero conc below read floor → Failed/insufficient_reads |
| `F_INSUF_P` | Only 3 non-zero concs pass floor → Failed/insufficient_points |
| `C0_ONLY` | Reads only at c=0 → Failed/insufficient_reads |

`F_LOW_R2` and `K_LOW` cover the lower edge of the clamp; `test_hill_pipeline.py` covers the upper one. Neither edge fails: an above-range row is clamped and flagged so the weak-binder tail stays in the table, and a below-range row gets its own `CensoredBelow` class so a block of equal bounds cannot head a downstream ranking.

The corpus runs at the shipped defaults — Hill exponent pinned at 1, hook detector off. The two cases whose intent needs a non-default knob (`P_N_HIGH`, `F_HOOK` / `NB_HOOK`) get their own opt-in tests in `test_corpus_e2e.py`.

## No-bin-mode clonotypes

| Key | Intent |
|---|---|
| `NB_GOOD` | Clean sigmoid frequency rise → Good |
| `NB_HOOK` | Top-1 freq < top-2/top-3 by > 0.02. Default run: detector off → Partial. With `hook_detector_enabled=True`: Failed/non_monotonic_signal |
| `NB_LOW_R2` | Scrambled frequencies → R² below `r2_threshold_failed`, reported Partial |
| `NB_FILLER` | Absorbs remaining sample depth so other clonotypes' frequencies reflect the intended profile. Not in manifest — tests ignore it. |

c=0 is deliberately absent from this corpus: `compute_global_baseline` takes the arithmetic mean of c=0 signals across clonotypes, and one clonotype dominating a c=0 sample would skew it for every other fit. Omitting c=0 forces the 4-parameter Hill fit, which is the production path for no-bin mode.

## Antigen-filter corpus

Four clonotypes, two labelled `antigen="target"` and two `antigen="other"`. Pipeline invoked with `target_antigen="target"` must drop distractors from all three output frames. Exercises R4.

## Regenerating

Run from `blocks/titeseq-analysis/software`:

```shell
uv run python tests/fixtures/generate_corpus.py
```

The generator is deterministic (`MASTER_SEED = 20260417`). Each clonotype uses a per-label derived sub-seed so adding a new entry doesn't ripple-shift earlier ones.

Commit the regenerated `*.parquet` files and `manifest.json` together with the code change that required the update.

## Manifest structure

Each clonotype entry has:

- `expected_class` or `expected_class_in` — single class or list of acceptable classes.
- `expected_reason` or `expected_reason_in` — single reason or list (for entries where the pipeline can land in more than one Failed branch without the test losing its intent).
- `kd_range` / `hill_range` — tolerance bands (not exact values) for Good/Partial fits. `null` for Failed clonotypes.
- `kd_out_of_range` — expected R14b flag, or `null` when not asserted.
- `kd_censor` — `"below"`, `"above"`, or `null`. Always asserted. No-bin mode never clamps, so every no-bin entry declares `null`.
- `r2_below` — present only on entries whose point is a bad fit that is reported anyway. Replaces the deleted `low_r2` verdict.
- `hill_plot_position_is_sentinel` — when `true`, R17 must emit `hillPlotPosition == -1.0` and `kdPlotPosition == max_conc * 10`.
- `notes` — human-readable reason this case exists.

Tolerance bands (not exact values) let scipy version drift without breaking the suite.
