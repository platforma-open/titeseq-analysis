# @platforma-open/platforma-open.titeseq-analysis

## 2.1.1

### Patch Changes

- Updated dependencies [f42932b]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.1.1
  - @platforma-open/platforma-open.titeseq-analysis.ui@2.2.1

## 2.1.0

### Minor Changes

- a3ad300: - Switch `params.json` encoding from the generic `json` module to
  `@platforma-sdk/workflow-tengo:canonical` so repeated runs produce
  byte-identical input. Eliminates a transient `CIDConflictError` path.
  - Bump `@platforma-sdk/*` catalog entries (workflow-tengo 5.15.0,
    block-tools 2.7.13, model/ui-vue/test 1.67.0, tengo-builder 2.5.14).
  - Expand `validationWarnings` to surface backend output errors
    (CID conflicts, exec failures) and mirror every `.args()` invariant as
    a severity:"error" issue — a disabled Run button now always has a
    visible reason.
  - Render validation alerts on every page (not just the settings drawer)
    via the shared `TiteseqPage` shell.
  - Enforce integer domains on `minReadsPerConcentration`,
    `minConcentrationPoints`, and `hookEffectMinReads` with inline
    "Must be a whole number" errors, and mark them `required` so the
    red-asterisk styling matches the ref-picker inputs.
  - Enable CSV export on the Clonotype Fit Results table via the
    SDK-native `show-export-button` prop, and add a new "Mean Bin Data"
    section that renders `signalPf` as a `PlAgDataTableV2` with the same
    export button. Together these cover the per-clonotype Hill fit
    outputs and the per-concentration mean-bin/fitted-mean-bin data the
    curves are derived from.
  - Realign the concentration axis to spec. Previously the workflow emitted
    a Long attomolar axis (`× 1e18`) that baked in a hidden molar
    assumption, violating R2 ("values are dimensionless floats") and
    rendering graph X-axis ticks at `10^6 … 10^12` (alien aM magnitudes).
    Now uses a String axis (canonical concentration string is the join
    key, per R14) plus a separate `concentrationValue:Double` sidecar
    PColumn that supplies the numeric source for log-scale graph
    rendering. Graph Maker plots `y = meanBin` against
    `x = concentrationValue`, joined on the shared String axis, so X-axis
    ticks render at the user's actual input concentrations. Drops the
    misleading `MAX_CONCENTRATION_M ≈ 9.2` rejection — any input unit /
    magnitude now works. Single import per TSV (no more dual-pcolumn
    `_Internal`/`_Export` variants), which also closes the dual-import
    CID-conflict pattern as a side effect. Spec calls for
    `concentration:Float` but the SDK gates axis types to
    `Int|Long|String` (see
    `core/platforma/sdk/workflow-tengo/src/pt/util.lib.tengo:352`).

### Patch Changes

- Updated dependencies [a3ad300]
  - @platforma-open/platforma-open.titeseq-analysis.workflow@2.1.0
  - @platforma-open/platforma-open.titeseq-analysis.model@2.1.0
  - @platforma-open/platforma-open.titeseq-analysis.ui@2.2.0

## 2.0.4

### Patch Changes

- @platforma-open/platforma-open.titeseq-analysis.workflow@2.0.1

## 2.0.3

### Patch Changes

- Updated dependencies [adacb76]
  - @platforma-open/platforma-open.titeseq-analysis.ui@2.1.1

## 2.0.2

### Patch Changes

- Updated dependencies [e9af4e2]
  - @platforma-open/platforma-open.titeseq-analysis.ui@2.1.0

## 2.0.1

### Patch Changes

- b374883: Update block description

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
  - @platforma-open/platforma-open.titeseq-analysis.ui@2.0.0
  - @platforma-open/platforma-open.titeseq-analysis.workflow@2.0.0
