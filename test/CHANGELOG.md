# @platforma-open/platforma-open.titeseq-analysis.test

## 2.0.0

### Major Changes

- d6d08b1: Initial public release of the Tite-Seq Analysis block.

  The block fits Hill curves per clonotype against MiXCR abundance across concentrations, emits Kd,app with a confidence class (Good / Failed), and renders Titration Curves, Kd Distribution, and Affinity vs Fit Quality tabs.

  **FACS sort-fraction correction.** An optional per-sample `sort_fraction` metadata column (C_bc/C_c from Adams, Mora, Walczak, Kinney 2016) sort-yield-corrects Mean Bin. When supplied, the output carries `pl7.app/titeseq/facsCorrected="true"`; absent the column, behaviour is bit-exact with the uncorrected pipeline. The Concentration and Sort-fraction dropdowns disambiguate by data — cross-exclusion plus a [0, 1]-range guard — so sort_fraction never appears under Antigen Concentration and vice versa.

  **Graph defaults.** Kd Distribution opens with a log y-axis and a green bin fill. Affinity vs Fit Quality opens with a log x-axis, filters to Failed clonotypes with Hill coefficient ≥ 0, and colours points by fit-failure reason (low_r2 → green, n_out_of_range → purple). Defaults bind to stable summaryPf columns so they transfer to any new block or project.

  **Spec compliance.** R5 (sample-metadata uniqueness), R9b (top-1/top-2 hook gate), R10 (Hill baseline bounds), R14 (dual concentration axes), R15 (Titration Curves default layout with affinity filter), R17 (Failed clonotypes park at Kd = -1.0 and render on Affinity vs Fit Quality).

  **R14 invariant.** `meanBin` and `fittedMeanBin` carry both `concentrationStr` (String join key) and `concentrationAM` (Long, attomolar, log-scale axis); a parametrized test guards the round-trip.

  **Inputs.** Empty, NaN, and Inf concentration columns are rejected; all-empty metadata columns are hidden from pickers; Antigen concentration accepts only Float or Double; the Target picker label follows the selected Antigen column.

  **UX.** The block subtitle derives from the first three populated inputs. Each tab shows a single page title. The Hill coefficient PColumn carries a Langmuir tooltip. Input labels and warnings are repositioned for clarity.

  **Fit Log.** Streams timestamped stage progress (load → validate → normalize → baseline → hook-detect → per-clonotype fit → write) and surfaces validator warnings.

  **Performance.** Hill fits parallelize above 50 clonotypes via `ProcessPoolExecutor`. Window functions replace self-joins in normalization. Validators run single-pass boundary counters. The metadata-uniqueness check scales with n_samples.

  **Robustness.** scipy `OptimizeWarning` marks fits Failed. The attomolar ceiling is guarded. Concentration casts defensively to Float64.

  **Architecture.** Migrated to `BlockModelV3`; concentration, bin, and antigen columns are anchored; Python runtime pinned to LTS-CPU wheels.

  **Tests.** Block-level integration exercises Samples-and-Data → import-vdj → Tite-Seq against a 28-sample synthetic fixture (5 clonotypes × 4 bins × 7 concentrations). A deterministic e2e corpus suite bounds numeric drift on per-clonotype columns.

### Patch Changes

- Updated dependencies [d6d08b1]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.0.0
