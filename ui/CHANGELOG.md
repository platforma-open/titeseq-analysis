# @platforma-open/platforma-open.titeseq-analysis.ui

## 2.5.2

### Patch Changes

- 3540e34: Keep the block subtitle in sync with the chosen input's current label
  - @platforma-open/platforma-open.titeseq-analysis.model@2.5.2

## 2.5.1

### Patch Changes

- 3e84ae6: Compute the block label without the UI

  The trace label is now the custom block label or "Tite-Seq Analysis (<target antigen>)",
  computed by the model, so opening the block no longer changes its arguments and a block
  that was never opened gets the same label. The abundance column's name moves to the
  subtitle only, written when the user picks the abundance; it was already part of the
  inherited upstream trace. Existing blocks show as stale once after the update because
  the arguments change shape.

- Updated dependencies [3e84ae6]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.5.1

## 2.5.0

### Minor Changes

- 43eecfb: Handle sample replicates in one block. Without a replicate column a block fits one titration:
  samples sharing a concentration pool into one point, and with the sort-fraction correction on the
  run stops with "sort_fraction values must sum to 1.0 per concentration", because two replicates sum
  to two.

  - New optional setting `replicateColumnRef` (block data and the kind's init params): a per-sample
    String or integer column naming each sample's replicate. Unset, every output and every byte of
    the args and the command line is the one the block produced before. One spec does change without
    a replicate column: the fit failure reason is now visible by default rather than hidden. A column
    that cannot be used — missing, empty, or two values on one sample — warns and runs as one
    replicate.
  - Each replicate is a titration of its own. The gate depth, the library depth, the sort-fraction
    sum and the missing (bin, concentration) check are taken within the replicate, so the replicates
    no longer divide each other's reads.
  - Each replicate is fitted alone; reads are never pooled across replicates. The combined value is
    the plain mean of log10 Kd,app. A replicate whose Kd,app ran past a wall of the titration enters
    pinned at that wall, as Adams's and Phillips's code both do.
  - New setting `replicateErrorModel`, shown as "Sample SD / √n" (the default) and
    "Population SD / √(n+1)". `sd_over_sqrt_n` is the sample SD over the square
    root of the replicate count, falling back to the single fit's own standard error for one
    replicate. `phillips` is the population SD over the square root of the count plus one, which is 0
    for one replicate; it reproduces the standard error published with Phillips 2021.
    `pl7.app/vdj/kdLog10Sd` reports the sample SD across the replicates beside it. Both are shown in
    the main table.
  - Two pre-combine filters, each a setting and each off unless set: `replicateFitErrorCutoff`
    (decades; Phillips 2021 published 1.0) and `replicateR2Floor` (Phillips 2021 published 0.8). Both
    apply to a replicate the fit measured inside the titrated range; a censored one is exempt. A
    filtered replicate is left out of the mean and counted, and a clonotype whose every replicate is
    filtered is Failed with reason `all_replicates_filtered`. The log fluorescence fit's R² gate is
    separate: it runs inside each replicate's own fit, and a gated replicate counts as failed.
  - `replicateCensoredMajority` (default 0.5): when no replicate measured a Kd,app, one censoring side
    must hold at least this share of the replicates and exceed the other side for that wall to be
    reported; otherwise the clonotype is Failed with reason `conflicting_censoring`. This is the
    block's own rule — no published method defines one. Under the log fluorescence fit the combined
    bound rates by the wall it sits at, `Below range` or `Above range`, as a single fit does.
  - The fit error cutoff reads a replicate's own fit standard error. The Mean bin fit reports none,
    so the setting is hidden and left out of the args under that equation, and a missing error no
    longer reads as a wide one. Sending it under the Mean bin fit previously rated every variant
    `Failed` with the reason `all_replicates_filtered`.
  - `pl7.app/vdj/replicateCount` and its family — `replicateMeasuredCount`,
    `replicateCensoredBelowCount`, `replicateCensoredAboveCount`, `replicateFailedCount`,
    `replicateFilteredCount` — count what each clonotype's replicates produced.
  - `pl7.app/vdj/kdLog10Se` is emitted under both equations when replicates are on: with several
    replicates the error is their spread, which does not depend on the equation.
  - A new frame on a (clonotype, replicate) axis, `replicatePf`, carries each replicate's own
    `pl7.app/vdj/kd`, `kdLog10Se`, `curveFitR2`, `affinityClass` and `kdReplicateState` (measured,
    censoredBelow, censoredAbove, failed, filtered). It is its own output, out of the exported bundle:
    the extra axis would fan out every row of the main table.
  - Clicking a variant in the main table opens its replicates beside it: each replicate's own
    Kd,app, standard error, R², class and fit failure reason. The panel's state and the replicate
    settings are new block data, added by a v3 -> v4 migration that moves no argument.
  - The signal frames gain the replicate on their key when a replicate column is set, so the
    Titration Curves graph draws each replicate's own curve instead of one pooled curve that belongs
    to no fit.

  A variant whose every replicate failed rates `Failed` with the reason `all_replicates_failed`, and
  each replicate's row in that panel carries its own fit failure reason. The count of
  replicates that entered the mean is labelled "Replicates in mean". The main table shows the fit
  failure reason by default.

  Under the "Population SD / √(n+1)" error model, a variant with one replicate in the mean reports that
  replicate's own fit error, not the paper's 0. A zero error rated every single-replicate variant
  Good and left the log-scale Affinity vs Fit Quality chart without a range.

- c6cb4c5: Choose the fit equation.

  - Two equations:
    - "Log fluorescence fit" (`log_fluorescence`), the new default: a Hill fit on log gate
      fluorescence with the exponent fixed at 1. Its one setting is the plateau floor, 0.667 by default. An optional
      column gives each gate's measured fluorescence. Without it, the gates are equally spaced.
      With 0 M samples, every clonotype shares one baseline, the read-weighted mean 0 M signal.
    - "Mean bin fit" (`hill_gate_index`): the 2.4.0 fit, unchanged.
  - Under `log_fluorescence` the block outputs `pl7.app/vdj/kdLog10Se`, the SE of log10 Kd,app. The
    Affinity vs Fit Quality page plots Kd,app against it, on a log y-axis. The Hill coefficient is a
    constant 1, hidden in the table and not exported. With measured gate fluorescence, the titration
    curves plot mean log10 fluorescence.
  - Each equation's workflow arguments carry only the settings that equation reads. A setting the
    drawer hides can neither disable Run nor make the block stale.
  - A migration keeps existing projects on `hill_gate_index`. Their workflow arguments do not change,
    so the upgrade does not make them stale.
  - The fit starts one worker per granted CPU, not one per host core.
  - Sort fractions are checked after the antigen filter, so a multi-antigen dataset passes the
    sum-to-1 check.
  - A clonotype with a constant signal is now Failed with no Kd,app. It was Partial, with a Kd,app set
    by the fit's amplitude bound.

- 3ea39c8: Replicate drill-down: on a run that combined sample replicates, double-clicking a variant in the
  main table opens a panel with that variant's replicates, one row each: Kd,app, SE, R², affinity
  class, fit failure reason, and the replicate state that says whether it entered the combined
  mean. Opening the panel does not stale the block or re-run anything.
- 7d02f2b: Log fluorescence fit:

  - A Kd,app that is a bound rates by the wall it sits at: `Below range` when the variant binds
    tighter than the lowest concentration titrated, `Above range` when it binds weaker than the
    highest. `Above range` covers a signal that stayed flat and dim and one that rose without
    reaching half-maximum. `Partial` now means an in-range Kd,app whose standard error is above 0.5.
    Lead Selection's default affinity cutoff keeps Good, Partial and Below range, so the library's
    tightest binders stay candidates and its weakest drop out.
  - "Kd,app out of range" carries the same fact without the side and moves to the optional columns.
  - An optional fit-quality gate on in-range fits: R² below its minimum rates `Failed`, `low_r2`,
    and, when a maximum is set, SE(log10 Kd,app) above it rates `Failed`, `high_kd_se`. Off by
    default. Phillips 2021 filtered at R² 0.8 and SE 1. Curve fit R² is hidden while the gate is
    off, and is empty instead of NaN for a signal with no variance.
  - An optional complete-titration rule: a clonotype with no reads at one of the titration's
    concentrations, 0 M included, rates `Failed`, `incomplete_titration`, before any fit. Off by
    default; Phillips 2021 fitted only complete titrations. Both fits read it.
  - "Min reads per concentration" and "Min concentration points" now apply to both fits, in one
    shared "Read coverage" section with the complete-titration rule. The log fluorescence fit used
    the same two values as fixed constants, so a block at the defaults fits exactly as before.
  - A new column, "-log10 Kd,app", reports Kd,app the way every published analysis of this assay
    does: a larger number is tighter binding, and the value shares a scale with the standard error
    beside it. It sits next to the affinity class, and Kd,app in molar moves to the optional columns.
  - Gate fluorescence above 7 log10 is refused as linear fluorescence. It used to move Kd,app by
    decades, or overflow the fit.
  - With no FACS bin column, the reason Run is disabled names the fix: choose the Mean bin fit.

  Both fits:

  - The main table shows one row per variant. The per-concentration signal no longer repeats each
    variant once per measured point; it stays on the Titration Curves page and in the export.
  - The missing (bin, concentration) warning counts the clonotypes of the selected antigen only, as a
    share of them, and no longer prints internal clonotype keys.
  - The Kd,app unit is the text in the concentration label's trailing parentheses: a column labelled
    "Antigen concentration (M)" gives "Kd,app (M)". Kd,app carries a description of what "apparent"
    means.

  The mean bin fit's columns and classes are unchanged.

  Fixes:

  - A project template now carries the gate fluorescence column, the maximum Kd,app SE and the
    complete-titration rule. A block created from a template that set them got none of the three.
  - An unsorted sample — reads at a concentration but no FACS bin — no longer drags the shared
    baseline below the first gate. Its reads were skipped in the mean's numerator while their weight
    stayed in the denominator, which pulled a 0 M baseline down by two thirds of a gate and moved
    flat non-binders into the titrated range as `Partial`. The mean bin fit keeps its own arithmetic,
    so its output is unchanged.

### Patch Changes

- db4e101: - Update graph-maker from 1.3.1 to 1.10.5. A saved chart now takes a default filter that lost an
  entry: a project last run under the mean bin fit no longer keeps the Failed-only filter on the
  Affinity vs Fit Quality page, which left that plot empty under the log fluorescence fit. The
  Kd Distribution tooltip no longer shows `Symbol(noGrouped)`.
  - Settings shows on the main page only. The other pages keep Logs.
- db4e101: Every number field in the settings drawer states each rule where the field itself checks it: its
  bounds are bounds, and a custom check covers only what a bound cannot say — a whole number, and the
  two paired limits (R² Good against Failed, Hill coefficient min against max). The bound messages
  come from the field component, so they read as they do in every other block.

  Run is disabled while the plateau floor is empty, as it already was for every other fit setting.

  A concentration column label with spaces no longer raises a warning. The label still becomes the
  Kd,app unit; spaces in it break nothing.

  Clearing a number field leaves it empty and marked "Value is required", and disables Run until a
  value is entered. A cleared box used to read empty over the number the run still used.

- e776c4a: UX review fixes:

  - The block overview describes both fit equations, the five affinity classes, replicate
    combination and the -log10 Kd,app column.
  - The "Fail in-range fits of poor quality" checkbox is gone. "Minimum R² (optional)" and "Maximum
    Kd,app SE (log10, optional)" always show under the log fluorescence fit, and each applies only
    when set. A new block starts at Minimum R² 0.8 and Maximum Kd,app SE 1. A project with the checkbox off
    keeps no limits; one with it on keeps both values.
  - "Require reads at every concentration" is no longer offered. Every run fits a titration with a
    gap.
  - "Min concentration points" accepts 5 or more. A stored project with 3 or 4 moves to 5, so it
    re-runs.
  - Kd,app in molar is shown again, right after the affinity class, in the main table and in the
    replicate panel.
  - Kd,app in molar shows three significant digits, for example 6.94e-9.
  - The mean bin fit also reports -log10 Kd,app, in the table, the replicate panel and the export.
  - Kd,app SE rounds to two decimals, like the SD. The per-replicate -log10 Kd,app rounds to two
    decimals, like the main table.
  - The SD column is "Kd,app SD (log10)". "Replicates measured" counts every replicate in the mean,
    including one above or below the titrated range. The other counts are "Replicates in range",
    "above range", "below range" and "failed". The old "Replicates measured" read 0 for a replicate at
    the top of the titration, which was measured.
  - The replicate panel title names the variant, and its Variant Id column is gone.
  - "Fit failure reason" is hidden by default in the main table. Columns shows it again.
  - The replicate panel no longer has "Entered the mean". It repeated the per-replicate affinity
    class.
  - Downstream blocks see four columns: affinity class, Kd,app, -log10 Kd,app and Kd,app SE. The
    replicate counts, the SD, R², the failure reason, the out-of-range flag and the per-concentration
    signal stay in this block's own table and graphs.
  - With replicates, double-clicking a variant row, or its Open button, opens that variant's
    replicates. The Open button no longer goes missing when the
    replicate outputs arrive after the table.
  - With the log fluorescence fit, "FACS bin" is marked required and says why Run is disabled while it
    is empty. "FACS sort fraction" and "Gate fluorescence" show before a bin is picked.
  - Log fluorescence fit: a Kd,app pinned at the edge of the range reports no R², as it reports no
    SE, and no bound shows an R² far below 0.
  - With replicates, the main table has no "Curve fit R²". Each replicate's R² shows in the replicate
    panel, as the published analyses report it.
  - "Equation" is now "Fit method equation", right above "FACS bin", since the equation decides
    whether a bin is required.
  - "Logs" shows on the Clonotype Fit Results page only, beside "Settings".
  - With replicates and a gate fluorescence column, each replicate is fitted against its own gate
    range. It was fitted against the range of all replicates together.

- Updated dependencies [43eecfb]
- Updated dependencies [c6cb4c5]
- Updated dependencies [db4e101]
- Updated dependencies [db4e101]
- Updated dependencies [3ea39c8]
- Updated dependencies [7d02f2b]
- Updated dependencies [e776c4a]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.5.0

## 2.4.0

### Minor Changes

- 03ca874: Update the SDK and accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) with the input read counts. The dataset list is built with buildDatasetOptions (requireEnrichments: false). The subset restricts the input before the fit, so the per-bin depths and clonotype frequencies are computed within the subset, and Kd,app values can differ from a full-data run. Outputs of a subset run carry pl7.app/inputSubset in their domain and a trace rooted in the subset, so their labels name it; the default subtitle names the subset too. Unfiltered runs are unchanged.

### Patch Changes

- Updated dependencies [03ca874]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.4.0

## 2.3.1

### Patch Changes

- e4569e0: Release every package through npm trusted publishing (OIDC) so each package's trusted publisher is used before it expires.
- Updated dependencies [e4569e0]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.2.1

## 2.3.0

### Minor Changes

- bccccbf: Accept synthetic-repertoire-profiler datasets as Tite-Seq input

  The read-count picker now also offers abundance anchors keyed on
  [sampleId, variantKey] (synthetic-repertoire-profiler's pl7.app/readCount),
  alongside the existing MiXCR [sampleId, clonotypeKey] / [sampleId, scClonotypeKey]
  anchors. The titration-curves plot picks its per-item key axis by excluding the
  concentration axis, so facet/grouping defaults resolve for variant datasets too.
  The workflow is unchanged — it was already axis-generic.

### Patch Changes

- Updated dependencies [bccccbf]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.2.0

## 2.2.3

### Patch Changes

- 3c3e66b: Fix labels machinery and CID conflicts for multi-instance TiteSeq blocks.

  - **Labels:** the workflow trace label now resolves to `customBlockLabel || defaultBlockLabel` (matching clonotype-clustering), so the user's subtitle override propagates to downstream pickers. UI auto-derivation leads with `targetAntigen` (the per-instance differentiator in multi-antigen studies) and drops the bin column — which was identical across instances and pushed the differentiator out of the join. Downstream blocks like Lead Selection now show distinct labels for each instance's "Affinity class" column instead of three identical entries.
  - **Determinism:** two TiteSeq blocks with byte-identical settings in the same project no longer hit `CIDConflictError`. The Python `fit-curves` binary now sorts every output frame before writing (`per_clonotype.tsv`, `fitted_mean_bin.tsv`, and a tie-breaker on `concentration_value.tsv`) and pins `maintain_order=True` on every `polars.group_by` whose result reaches a disk write. A prior commit handled `mean_bin.tsv` only; the same hazard lived in the three sibling outputs.

## 2.2.2

### Patch Changes

- Updated dependencies [d2ecf8c]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.1.2

## 2.2.1

### Patch Changes

- f42932b: fix additional curves
- Updated dependencies [f42932b]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.1.1

## 2.2.0

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
  - @platforma-open/platforma-open.titeseq-analysis.model@2.1.0

## 2.1.1

### Patch Changes

- adacb76: Remove whitespace in graph pages

## 2.1.0

### Minor Changes

- e9af4e2: Add fitted mean bin overlay to the Titration Curves page. When the Hill fit is available, the fitted curve renders alongside observed per-clonotype dots, faceted by clonotype. Uses the new `additionalCurves` slot introduced in `@milaboratories/graph-maker` 1.3.0.

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
