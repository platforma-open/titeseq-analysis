# @platforma-open/platforma-open.titeseq-analysis.kind

## 1.2.0

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

## 1.1.0

### Minor Changes

- 03ca874: Update the SDK and accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) with the input read counts. The dataset list is built with buildDatasetOptions (requireEnrichments: false). The subset restricts the input before the fit, so the per-bin depths and clonotype frequencies are computed within the subset, and Kd,app values can differ from a full-data run. Outputs of a subset run carry pl7.app/inputSubset in their domain and a trace rooted in the subset, so their labels name it; the default subtitle names the subset too. Unfiltered runs are unchanged.
