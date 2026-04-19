# Titeseq Analysis — Block Test Plan

**Scope.** All plausible block-level tests (`test/src/wf.test.ts`) to verify the titeseq-analysis block against the [affinity-profiling specs](../../../docs/text/work/projects/affinity-profiling/README.md) (R1–R19b) and the [PColumn spec](../../../docs/text/work/projects/affinity-profiling/pcolumn-spec.md). This is a menu of possibilities — not all will be implemented. Priorities (P0 must-have, P1 should-have, P2 nice-to-have) are provided to aid triage.

**Test categories.**
- **BT** — Block integration test via `blockTest` from `@platforma-sdk/test` (real backend, full wf + model + UI).
- **TPL** — Workflow template test via `tplTest` (Tengo in isolation, no full backend).
- **MDL** — Model unit test (TypeScript only: validation, argsValid, output derivation).
- **PY** — Python software test (in `software/tests/` — already present; referenced only).

Focus below is on **BT**, since `test/src/wf.test.ts` is currently empty.

---

## Review notes (2026-04-19)

Coverage of R1–R19b looks comprehensive. Gaps and clarifications worth addressing before implementation starts:

- **MDL harness is unspecified.** `@platforma-sdk/test` ships `blockTest` and `tplTest` but no MDL-style harness. Tests tagged **MDL** (T1.3, T5.1, T5.2, T6.1–T6.3, T11.1–T11.4) need either a decision to (a) run them via plain vitest with a hand-rolled mock `ctx`, or (b) roll them into **BT** tests against a running backend. Pick one before picking these up — otherwise each implementer will choose differently.
- **I1(b) is the critical path.** Every test beyond T1.1 depends on a synthetic upstream PColumn helper. Without it, the MVP collapses to one test. Treat I1(b) as its own ticket landing before T2.x.
- **Model export mismatch with I3 boilerplate.** `model/src/index.ts` exports `model`; the I3 snippet (and `blocks/mixcr-amplicon-alignment`) expects `platforma`. Add `export { model as platforma }` to the model so the canonical `import { ..., platforma } from 'this-block.model'` pattern works.
- **T1.1 is under-specified.** Only checks the four option arrays. See the tightened checklist below — also assert `validationWarnings: []`, `binMode: false`, and no workflow-run side effects.
- **T3.0 is missing.** A P0 "all three output handles present" smoke should precede T3.1+. It catches workflow return-shape regressions (e.g., missing `exports` key) faster than column-level assertions.
- **`validationWarnings` reactivity isn't covered.** T6.1–T6.3 assert a warning appears but not that it clears when the user fixes the input. Clearing behavior is where ordering bugs bite. Added as T6.1b–T6.3b.
- **`binMode` output isn't tested.** It drives R19b label swaps. Added as T11.5.
- **Log-handle assertion pattern isn't defined.** T6.4 says "surfaced via logHandle" — I3 should include a helper that reads `logHandle` content (via `ml.driverKit.blobDriver.getContent`) so every error-path test doesn't reinvent the wheel.

---

## 0. Test Infrastructure Prerequisites

Before meaningful block tests exist, the harness needs:

- [ ] **I1.** Upstream fixture provider. Block tests need an abundance PColumn with axes `[sampleId][clonotypeKey]` plus per-sample metadata (concentration, bin, antigen). Two options:
  - **(a) Live upstream chain** — `Samples & Data` → `MiXCR Clonotyping` (or synthetic upstream) wired into `setupProject()` as in `blocks/clonotype-browser/test/src/wf.test.ts` and `blocks/mixcr-amplicon-alignment/test/src/wf.test.ts`. Highest fidelity; slowest.
  - **(b) Synthetic PColumn helper** — a utility that materialises an abundance PColumn and per-sample numeric/int/string columns directly into the result pool without running MiXCR. Fast and deterministic; decouples titeseq-analysis tests from MiXCR versioning. **Recommended.**
- [ ] **I2.** Corpus reuse. The Python corpus in `software/tests/data/corpus/` (`reads_bin_mode.parquet`, `reads_no_bin_mode.parquet`, `reads_antigen.parquet`, `manifest.json`) already encodes the three canonical input scenarios. A block-test helper should be able to load these and emit them as a synthetic PColumn to avoid regenerating ground truth.
- [ ] **I3.** `blockTest` boilerplate — import `BlockArgs`, `BlockOutputs`, `platforma` from this-block; `awaitStableState`, `blockTest` from `@platforma-sdk/test`; `wrapOutputs`, `InferBlockState` from `@platforma-sdk/model`. Also add a `readLogHandle(ml, handle)` helper that wraps `ml.driverKit.blobDriver.getContent` so T6.4 and any other log-asserting test has one place to update.
- [ ] **I4.** MDL harness decision. `@platforma-sdk/test` exposes no MDL-style entry point. Decide between (a) plain vitest + hand-rolled mock `ctx` (fastest, but skips the real result-pool queries), or (b) promote every MDL test to BT so it runs against the live model. Write it up in this file and pick one approach for the whole suite.
- [ ] **I5.** Model exports `platforma` alias. Add `export { model as platforma }` to `model/src/index.ts` so the I3 boilerplate matches the `blocks/mixcr-amplicon-alignment` convention.

---

## 1. Tier 1 — Smoke and discovery (P0)

These exercise the model-layer option discovery and "empty project" path without running the Python fitter. Fast, cheap, and catch 80% of refactor regressions.

- [ ] **T1.1 (BT, P0).** *Empty project, no inputs.* Add the block to an empty project. Expect:
  - `abundanceOptions`, `concentrationOptions`, `binOptions`, `antigenOptions` each resolve to `{ ok: true, value: [] }`.
  - `validationWarnings` resolves to `{ ok: true, value: [] }` (defaults satisfy all validators).
  - `binMode` resolves to `{ ok: true, value: false }` (no `binColumnRef`).
  - `summaryPf`, `signalPf`, `logHandle` are absent or undefined — the `argsValid` gate should block the workflow from running.
  - `isRunning` is not truthy.

  Pattern: mirror `blocks/mixcr-amplicon-alignment/test/src/wf.test.ts`'s `empty inputs` test.
- [ ] **T1.2 (BT, P0).** *Options surface after upstream connects.* Add upstream fixture (I1b). Expect `abundanceOptions` to include the synthetic PColumn, `concentrationOptions` to include the numeric sample-metadata column, `binOptions` to include only integer sample-metadata, `antigenOptions` to include only string sample-metadata. Validates model output selectors in `model/src/index.ts` lines ~78–140.
- [ ] **T1.3 (MDL, P0).** *`argsValid` truth table.* Directly call the model with crafted `ctx`:
  - `abundanceRef` unset → invalid.
  - `concentrationColumnRef` unset → invalid.
  - Both set, no antigen column → valid.
  - `antigenColumnRef` set, `targetAntigen` empty → invalid (R4).
  - All required set with bin + antigen + targetAntigen → valid.

## 2. Tier 2 — Core happy paths (P0)

One test per canonical mode. These are the backbone of the suite.

- [ ] **T2.1 (BT, P0).** *Bin mode, 0 M control, single antigen.* Inputs: 4 bins, 6 log-spaced concentrations including 0 M, single antigen (no antigen column). Assert:
  - Workflow reaches stable state without error.
  - `summaryPf` contains 6 PColumns: `pl7.app/vdj/kd`, `hillCoefficient`, `curveFitR2`, `affinityClass`, `fitFailureReason`, `kdOutOfRange` — each with axis `[clonotypeKey]`, domain `{"pl7.app/block": blockId}`.
  - `signalPf` contains `meanBin` + `fittedMeanBin` with axes `[clonotypeKey][concentration:Float]`.
  - `exportPf` (pf) is populated (downstream Lead Selection consumer).
  - At least one clonotype classified `Good`; K_D,app within tolerance of ground truth per Milestone 1 (≥90% within 10%).
  - `pl7.app/unit` annotation on `kd` equals the concentration column's label (R2).
- [ ] **T2.2 (BT, P0).** *No-bin mode (R7b).* Same dataset minus bin column. Assert:
  - meanBin label swaps to `"Clonotype frequency"`, fittedMeanBin to `"Fitted clonotype frequency"` (R19b, see `main.tpl.tengo` lines labelling).
  - Signal values are frequencies (0–1 range) not bin indices.
  - Fits still produce K_D,app values (different hook-effect thresholds apply).
- [ ] **T2.3 (BT, P0).** *Antigen column provided.* Multi-antigen dataset (say 2 antigens × 4 bins × 6 concentrations). Set `targetAntigen` to antigen A. Assert:
  - Only antigen-A samples reach the fitter (R4 filter, Stage 1 of processing pipeline).
  - Fits match T2.1 for antigen-A samples.
  - Switching `targetAntigen` to antigen B without rerunning upstream yields a different fit set.
- [ ] **T2.4 (BT, P0).** *No 0 M control → 4-parameter fit.* Same as T2.1 without the 0 M concentration point. Assert successful fits and that `kd` values for the same clonotypes differ from T2.1 (4-param vs 3-param baseline) but remain within spec tolerance.

## 3. Tier 3 — Output-shape invariants (P0–P1)

Asserts on the structural contract between the block and downstream consumers. Cheap to check alongside T2.x.

- [ ] **T3.0 (BT, P0).** *All three output handles present.* After a minimal successful run (can piggyback on T2.1), assert that `stableState.outputs` resolves `summaryPf`, `signalPf`, and `logHandle` with `ok: true`. Catches workflow-shape regressions — e.g., a missing `exports` key or an unreturned output — before any column-level assertion.
- [ ] **T3.1 (BT, P0).** *PColumn axes and value types match `pcolumn-spec.md`.* For each summary column, assert:
  - `kd`, `hillCoefficient`, `curveFitR2` → `valueType: Double`, axes `[clonotypeKey]`.
  - `affinityClass`, `fitFailureReason` → `valueType: String`, axes `[clonotypeKey]`.
  - `kdOutOfRange` → `valueType: String` (boolean-as-string), axes `[clonotypeKey]`.
  - `meanBin`, `fittedMeanBin` → `valueType: Double`, axes `[clonotypeKey][concentration:Float]`.
- [ ] **T3.2 (BT, P0).** *Key annotations present on each output.* Assert `pl7.app/label`, `pl7.app/isScore`, `pl7.app/table/visibility`, `pl7.app/table/orderPriority` match the values in `pcolumn-spec.md`. Critically: `affinityClass` carries `pl7.app/isDiscreteFilter: "true"` and `pl7.app/discreteValues: '["Good","Partial","Failed"]'` (Lead Selection integration).
- [ ] **T3.3 (BT, P0).** *`pl7.app/unit` propagates from concentration column label.* Set concentration column label to `"nM"`; assert `kd.annotations["pl7.app/unit"] === "nM"`. Change to `"µM"`; assert it updates. Omit label; assert the annotation is absent.
- [ ] **T3.4 (BT, P0).** *c=0 excluded from meanBin/fittedMeanBin output (R14).* Include 0 M in inputs; assert the concentration axis of both signal PColumns contains no `0.0` key (log-scale rendering would break otherwise).
- [ ] **T3.5 (BT, P1).** *`fittedMeanBin` is null for Failed-fit clonotypes.* Inject a deliberately un-fittable clonotype (all zeros after floor); assert its `fittedMeanBin` rows are missing/null, but its `meanBin` rows are present (so users see the raw dots without a curve).
- [ ] **T3.6 (BT, P1).** *`pl7.app/trace` provenance.* Assert trace annotations on output columns chain back to the source abundance column (via `pSpec.makeTrace` in `main.tpl.tengo`).
- [ ] **T3.7 (BT, P1).** *`blockId` domain scoping.* All outputs carry `domain.pl7.app/block == blockId` so two instances of the block don't collide in the result pool.

## 4. Tier 4 — Classification matrix (P1)

Targets `classify.py` truth table via black-box assertions on the block output. Each input is crafted to land on a specific row.

- [ ] **T4.1 (BT, P1).** *Good class.* Clean synthetic curve with R²_w ≥ 0.8 and n ∈ [0.5, 2.0] → `affinityClass == "Good"`, `fitFailureReason == null`, `kd` populated.
- [ ] **T4.2 (BT, P1).** *Partial (medium R²).* Moderate noise so 0.5 ≤ R²_w < 0.8, n in range → `Partial`, `fitFailureReason == null`.
- [ ] **T4.3 (BT, P1).** *Partial (Good R² downgraded by n).* R²_w ≥ 0.8 but n < `nMin` (e.g., 0.3) → `Partial`, `fitFailureReason == "n_out_of_range"`.
- [ ] **T4.4 (BT, P1).** *Failed (low R²).* Pure noise → R²_w < 0.5 → `Failed`, `fitFailureReason == "low_r2"`.
- [ ] **T4.5 (BT, P1).** *Failed (Partial R² with bad n).* 0.5 ≤ R²_w < 0.8 and n out of range → `Failed`, `fitFailureReason == "n_out_of_range"`.
- [ ] **T4.6 (BT, P1).** *Failed (non-monotonic / hook effect, R9b).* Clonotype whose mean bin drops at the max concentration by > `hookEffectThresholdBin` → `Failed`, `fitFailureReason == "non_monotonic_signal"`, `fittedMeanBin` null, `kd` null.
- [ ] **T4.7 (BT, P1).** *Failed (insufficient_reads).* All concentration points below `minReadsPerConcentration` → `Failed`, `fitFailureReason == "insufficient_reads"`.
- [ ] **T4.8 (BT, P1).** *Failed (insufficient_points).* Reads survive at exactly `minConcentrationPoints − 1` points → `Failed`, `fitFailureReason == "insufficient_points"`.
- [ ] **T4.9 (BT, P2).** *Failed (convergence_failure).* Construct a curve scipy cannot fit (e.g., all points at the same bin value giving amplitude ≈ 0, or adversarially noisy) → `Failed`, `fitFailureReason == "convergence_failure"`.
- [ ] **T4.10 (BT, P1).** *Mixed population.* Single test with a curated mix (20 clonotypes spanning all 6 failure reasons + Good + Partial). Assert the exact count per class — single-assert proxy for the full matrix if T4.1–T4.9 are too slow.

## 5. Tier 5 — Parameter-change reactivity (P1)

R19 mandates that viz parameters re-derive without rerunning the workflow. These tests verify that the right knobs live in the model layer vs the Python layer.

- [ ] **T5.1 (MDL, P1).** *Changing `r2ThresholdGood` downstream of a run re-partitions classes.* Run once, capture outputs, change threshold, assert classifications flip without workflow rerun (requires classification to be a model-derived output — current implementation re-runs Python; if so, mark this test as a design-question probe).
- [ ] **T5.2 (MDL, P1).** *Changing `nMin`/`nMax` downstream changes class.* Same pattern.
- [ ] **T5.3 (BT, P1).** *Changing `minReadsPerConcentration` forces rerun.* Floor filter is upstream of signal computation → workflow must rerun. Assert `isRunning` transitions true→false on arg change.
- [ ] **T5.4 (BT, P2).** *Changing `customBlockLabel` does not trigger a rerun.* Label-only state; assert workflow doesn't re-execute (catch spurious reactivity).

## 6. Tier 6 — Validation and user errors (P1)

Cover the `validationWarnings` path in `model/src/index.ts` and upstream workflow-level errors.

- [ ] **T6.1 (MDL, P1).** *Concentration column with spaces warns (R2).* Label `"Antigen concentration (nM)"` → warning in `validationWarnings`.
- [ ] **T6.1b (MDL, P1).** *Warning clears when label is fixed.* After T6.1 fires, switch the concentration-column label to `"nM"`; assert the warning disappears from `validationWarnings` (reactivity).
- [ ] **T6.2 (MDL, P1).** *Antigen column set but `targetAntigen` empty → error severity.*
- [ ] **T6.2b (MDL, P1).** *Error clears once `targetAntigen` is set.* Assert `argsValid` flips to true and the error disappears once `targetAntigen` is non-empty.
- [ ] **T6.3 (MDL, P1).** *`targetAntigen` set but no antigen column → warning (value ignored).*
- [ ] **T6.3b (MDL, P1).** *Warning clears when `targetAntigen` is cleared.* Assert the "value ignored" warning disappears once `targetAntigen` is reset to `""`.
- [ ] **T6.4 (BT, P1).** *Negative concentration rejected (R2).* Inject a metadata column with a `-1` value; expect a workflow-level error surfaced via `logHandle` or a structured error output, not a silent bad fit.
- [ ] **T6.5 (BT, P2).** *Zero concentration without bin column produces a warning (R2).* Ambiguous 0 M control in no-bin mode.
- [ ] **T6.6 (BT, P1).** *Non-integer bin column rejected (R3).* If `binColumnRef` points to a Float column, the model option filter should exclude it; assert it never appears in `binOptions`.

## 7. Tier 7 — Edge cases by data shape (P1–P2)

- [ ] **T7.1 (BT, P1).** *Minimum viable dataset.* Exactly `minConcentrationPoints` + `minReadsPerConcentration` everywhere → all fits succeed, no failures for threshold reasons.
- [ ] **T7.2 (BT, P1).** *One sample per concentration.* Bin mode with only one bin sample per concentration — degenerate but not illegal. Should produce a result (possibly Failed via `insufficient_points`).
- [ ] **T7.3 (BT, P2).** *8-bin FACS dataset.* R3 explicitly notes 8-bin is common. Verify the pipeline doesn't hard-code 4 bins anywhere — run with bin labels `1..8`.
- [ ] **T7.4 (BT, P2).** *Bin labels starting from 0.* Python `max_bin_label` and weight normalization should handle `{0,1,2,3}` as well as `{1,2,3,4}`.
- [ ] **T7.5 (BT, P2).** *Very large clonotype count (scalability).* 10k+ clonotypes — exercises `processColumn`/`ProcessPoolExecutor`, memory, and `tsvFileBuilder` at scale. May be a regression-only test, not in CI.
- [ ] **T7.6 (BT, P2).** *Single clonotype.* Degenerate smallest case. Ensures no div-by-N-clonotypes anywhere.
- [ ] **T7.7 (BT, P1).** *Single concentration point.* Should classify all clonotypes as `insufficient_points`; no crash.
- [ ] **T7.8 (BT, P1).** *Highly unbalanced read coverage.* One concentration point dominates, others near floor. Tests `σ_j = 1/√w_j` weighting gives reasonable results (vs uniform fit).
- [ ] **T7.9 (BT, P2).** *Duplicate sample IDs in metadata.* R5 validation path — assert `validate_sample_metadata_uniqueness` surfaces a user-readable error.
- [ ] **T7.10 (BT, P2).** *Missing metadata for some samples.* Concentration column doesn't cover all samples in the abundance column. Assert either rejected at validation or only matching samples processed.

## 8. Tier 8 — Concentration axis canonicalization (P0)

R14's "concentration axis key canonicalization" is a high-risk regression area — float-serialization drift between TS → Tengo → Python → Parquet is real.

- [ ] **T8.1 (BT, P0).** *Concentrations with many decimals round-trip intact.* Use `0.001`, `1e-9`, `0.0000001`; assert that `meanBin` and `fittedMeanBin` share identical concentration axis keys (can be joined in Graph Maker without drift).
- [ ] **T8.2 (BT, P1).** *Scientific-notation string stability.* Metadata label `1.0e-6` vs `0.000001` should yield a consistent axis key; document which canonical form wins.
- [ ] **T8.3 (BT, P1).** *Join integrity between meanBin and fittedMeanBin.* Build a test that consumes both PColumns via the Graph Maker helper and asserts every `(clonotypeKey, concentration)` point in `fittedMeanBin` has a matching point in `meanBin`. This is the real user-facing failure mode if canonicalization breaks.

## 9. Tier 9 — `kdOutOfRange` flag (P1)

R14b requires K_D values outside the experimental concentration range to be flagged but **not** suppressed.

- [ ] **T9.1 (BT, P1).** *K_D below min concentration.* Tight binder saturating at the lowest point → `kd` value present, `kdOutOfRange == "true"`, `affinityClass` honored (Good/Partial/Failed independently).
- [ ] **T9.2 (BT, P1).** *K_D above max concentration.* Weak binder not saturating → `kd` extrapolated, `kdOutOfRange == "true"`.
- [ ] **T9.3 (BT, P1).** *K_D inside range.* Normal case → `kdOutOfRange == "false"`.

## 10. Tier 10 — Numerical regression (P1)

Anchors the block against drift in scipy / polars / numpy / Python version upgrades. The software-side `tests/regression/test_synthetic_titeseq.py` handles deep fitting regression; block tests add end-to-end coverage.

- [ ] **T10.1 (BT, P1).** *Golden-output regression.* For the canonical bin-mode corpus (reused from `software/tests/data/corpus/reads_bin_mode.parquet`), snapshot `kd`, `hillCoefficient`, `r2`, `affinityClass` per clonotype. Assert within tight numerical tolerance (e.g., 1e-4 relative) on subsequent runs. Protects both the Python fitter and the Tengo→Parquet→PColumn import path.
- [ ] **T10.2 (BT, P2).** *Deterministic under parallelism.* Run with `--workers 1`, `--workers 4`, `--workers auto`; assert outputs are byte-identical (or within tolerance if non-associative float ops are unavoidable).
- [ ] **T10.3 (BT, P1).** *Recovery of known K_D values (Milestone 1 acceptance).* Generate 50 clonotypes with known K_D spanning 0.1×–100× of the concentration grid; assert ≥90% recover within 10% of truth. This is the Milestone 1 test called out verbatim in the spec.

## 11. Tier 11 — UI / model-output smokes (P2)

Quick sanity that model-side derivations the UI consumes actually populate. These are cheap proxies for running the UI itself.

- [ ] **T11.1 (MDL, P2).** *`displayTable` output has rows.* After a successful run, `ctx.outputs.resolve('tableData')` returns a non-empty handle; `createPlDataTableV2` does not throw.
- [ ] **T11.2 (MDL, P2).** *Graph Maker PFrames build.* The titration-curves and KD-histogram `PFrameHandle` outputs (via `createPFrameForGraphs`) are non-null and contain the expected columns.
- [ ] **T11.3 (MDL, P2).** *`isRunning` transitions.* Block state shows `isRunning: true` during execution, `false` once outputs stabilize.
- [ ] **T11.4 (MDL, P2).** *`customBlockLabel` in title.* Set a custom label; assert `ctx.title()` reflects it (R19-adjacent UX).
- [ ] **T11.5 (MDL, P1).** *`binMode` output flips with `binColumnRef`.* Assert `binMode` is `false` when `binColumnRef` is undefined and `true` when set. Drives R19b label swaps in the UI — regressions here silently break the "Mean bin" vs "Clonotype frequency" rendering.

## 12. Tier 12 — Cross-block integration (P2)

- [ ] **T12.1 (BT, P2).** *Lead Selection consumption.* Chain an `antibody-tcr-lead-selection` block downstream (if available as a workspace dep); assert the In Vitro preset picks up `kd` and `affinityClass` automatically. Validates R13 + spec's stated downstream connection.
- [ ] **T12.2 (BT, P2).** *Re-running upstream changes outputs.* Modify the upstream abundance PColumn (e.g., add a sample); assert titeseq re-runs and outputs update.
- [ ] **T12.3 (BT, P2).** *Multiple titeseq blocks in one project.* Two instances targeting different antigens from the same upstream — `blockId` domain scoping should prevent PColumn collisions (covered by T3.7 but worth as an end-to-end case).

## 13. Tier 13 — Workflow template isolation (P2)

Cheap-ish `tplTest` cases that don't need a full backend. Useful when iterating on Tengo logic.

- [ ] **T13.1 (TPL, P2).** *Main template produces three expected output file handles.* `per_clonotype.tsv`, `mean_bin.tsv`, `fitted_mean_bin.tsv` are declared as `saveFile` outputs in the exec builder (`main.tpl.tengo` Stage E).
- [ ] **T13.2 (TPL, P2).** *`setTableProps` helper applies annotations correctly.* Pure Tengo unit.
- [ ] **T13.3 (TPL, P2).** *`hasBin`/`hasAntigen` branching.* Template should omit the bin header and antigen column when refs are undefined.

## 14. Tier 14 — Performance / resource (P2)

Not part of CI. Run manually when touching memory-sensitive stages.

- [ ] **T14.1 (BT, P2).** *16 GiB TSV builder limit.* Generate an abundance PColumn close to the `readsBuilder.mem("16GiB")` bound in Stage B and verify graceful handling.
- [ ] **T14.2 (BT, P2).** *`cpu: "auto"` parallelism scaling.* Measure Python stage wall-time at different `--workers` values.

---

## Recommended MVP (what to implement first)

If picking a minimum viable set to ship alongside Milestone 1, I'd take:

1. **T1.1, T1.2** — zero-cost empty/options smokes. *T1.1 is the only MVP test implementable before I1(b) lands.*
2. **T2.1, T2.2, T2.3, T2.4** — the four canonical mode combinations.
3. **T3.0, T3.1, T3.2, T3.3, T3.4** — output contract for downstream consumers.
4. **T4.10** — mixed-population classification matrix (subsumes T4.1–T4.9).
5. **T8.1, T8.3** — axis canonicalization (high regression risk).
6. **T10.3** — Milestone 1 acceptance test (K_D recovery within 10%).

That's ~13 tests, covers every requirement R1–R14b end-to-end, and gives a regression-safety net for R15–R19 via output-shape assertions even before the Graph Maker overlay slot lands.

---

## Cross-reference: requirements coverage

| Req   | Covered by |
|-------|------------|
| R1    | T1.2, T2.1 |
| R2    | T3.3, T6.1, T6.4, T6.5 |
| R3    | T2.1, T6.6, T7.3, T7.4 |
| R4    | T2.3, T6.2, T6.3 |
| R5    | T6.4, T7.9, T7.10 |
| R6    | T2.1 vs T2.4 (delta) |
| R7    | T2.1, T8.3 |
| R7b   | T2.2 |
| R8    | T7.8 |
| R9    | T4.7, T4.8, T7.1 |
| R9b   | T4.6 |
| R10   | T10.3 |
| R11   | T4.1–T4.5 (R² thresholds) |
| R12   | T4.1–T4.10 |
| R13   | T3.1, T3.7, T12.1 |
| R14   | T3.4, T3.5, T8.1, T8.2, T8.3 |
| R14b  | T9.1, T9.2, T9.3 |
| R15   | T11.2 (PFrame build smoke); full viz tests pending Graph Maker slot |
| R16   | T11.2 |
| R17   | T11.2 |
| R18   | T11.1 |
| R19   | T5.1, T5.2, T5.4 |
| R19b  | T2.2 |

---

## Out of scope for this plan

- Python-internal unit tests (already in `software/tests/unit/`).
- UI component rendering tests (no Vue Testing Library harness currently wired).
- Graph Maker curve overlay slot behavior — blocked on Milestone 2 (visualization-api).
- Absolute K_D calibration, cross-experiment comparison — out of spec scope.
