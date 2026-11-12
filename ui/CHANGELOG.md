# @platforma-open/platforma-open.titeseq-analysis.ui

## 3.0.1

### Patch Changes

- @platforma-open/platforma-open.titeseq-analysis.model@3.0.1

## 3.0.0

### Major Changes

- 956a634: Publish the fit software again under a new version. The 2.3.0 archive in the software registry
  holds an older build, because that name was taken before the release and the upload was skipped.
  The 2.6.0 workflow passes options that build does not know, so a run on a backend that uses the
  archive stops with "unrecognized arguments".

### Patch Changes

- Updated dependencies [956a634]
  - @platforma-open/platforma-open.titeseq-analysis.model@3.0.0

## 2.6.0

### Minor Changes

- 68af0de: Report one Kd,app per amino-acid variant on a nucleotide variant dataset, with an error from its
  synonymous nucleotide variants, and fix the error across sample replicates.

  Synonymous variants:

  - Turns on by itself when the dataset's variant axis is nucleotide and its nt -> aa linker
    (`pl7.app/repertoire/aaToNt`) resolves. Every other input runs as in 2.5.0.
  - Each nucleotide variant is fitted on its own, exactly as in 2.5.0. Reads of two variants are never
    pooled. The fits are then combined per amino-acid variant: the mean of their log10 Kd,app. With
    two or more synonymous variants, the SE (`pl7.app/vdj/kdLog10Se`) is the larger of their spread,
    SD / sqrt(n), and their fit error, RMS(se_i) / sqrt(n), where the equation reports a fit SE. n
    counts every variant in the mean, so one with no fit SE, such as one at an edge, counts at the
    typical error of the others and cannot lower the SE. With one variant, its own fit SE where the
    equation reports one. A variant
    censored at the edge of the titration enters at that edge; a Failed one is left out and
    counted. With no variant measured inside the range, the value is the edge that holds at least
    half of them and more than the other edge, as replicates are combined; otherwise the amino-acid
    variant is Failed (`conflicting_censoring`), since a mean of opposite edges would land inside the
    range.
  - The amino-acid grain has a "Fit failure reason" column (`pl7.app/vdj/fitFailureReason`, optional
    in the table, exported): `all_variants_failed` when no synonymous variant has a Kd,app,
    `conflicting_censoring`, or `not_productive` for a variant that is not productive. With sample
    replicates it reads `all_replicates_failed`, `conflicting_censoring` or `not_productive`, and each
    amino-acid row of the per-replicate table carries its own reason.
  - The mean is always the plain mean, as Starr combines synonymous variants; there is no
    read-weighted option. Read-weighting let one deep synonymous variant swamp the others: two
    variants 2 decades apart, with 1,000 and 10 reads and no fit SE, rated Good with an SE of 0.04.
  - `pl7.app/vdj/kdLog10SynonymSdAll` on the amino-acid axis ("Kd,app SD across synonymous
    variants, edges included (log10)"), without replicates: the SD of log10 Kd,app across the
    synonymous variants in the mean, those at an edge entering at the edge, with no fit-error floor.
    Empty with fewer than two variants, for a bound and for a stop or failed variant. It has a name
    of its own because `pl7.app/vdj/kdLog10Sd` is the spread across replicates on every grain.
  - `pl7.app/vdj/kdLog10SynonymSd` on the amino-acid axis ("Kd,app SD across synonymous variants,
    measured in range only (log10)", optional in the table): the SD of log10 Kd,app across the
    synonymous variants measured inside the range; with replicates, each combined over the
    replicates first. With replicates it stays, beside the replicates' `kdLog10Sd`.
  - Amino-acid columns sit on the linker's amino-acid axis and reuse the nucleotide names:
    `pl7.app/vdj/kd`, `pl7.app/vdj/kdNegLog10`, `pl7.app/vdj/kdLog10Se`, `pl7.app/vdj/affinityClass`,
    `pl7.app/vdj/kdOutOfRange`. The affinity class takes the equation's declared values and cutoff,
    as on the nucleotide grain. `pl7.app/vdj/kdNtVariants` counts the nucleotide variants behind
    each value. The total carries the Kd,app columns' domain; the counts by state add
    `pl7.app/vdj/kdNtVariantState` (measured, censored low, censored high, failed). With replicates
    the total counts a variant fitted in several replicates once, the counts by state are summed
    over the replicates, and their descriptions say so.
  - `pl7.app/vdj/variantClass` (String, discrete filter) on both grains, after Starr's classes:
    `wildtype`, `synonymous`, `nonproductive`, `1 nonsynonymous`, `>1 nonsynonymous`. The amino-acid
    grain needs the amino-acid `pl7.app/repertoire/mutationCount`; the nucleotide grain needs that and the
    nucleotide one. `nonproductive` needs the amino-acid `pl7.app/vdj/sequence/productive`; without it there
    is no such class. Each input is used only when exactly one column matches it.
  - An amino-acid variant that is not productive (a stop codon or a frame shift) has no Kd,app, SE
    or parent delta and is rated Failed; the run
    summary, the Fit log and the block's warnings count them. Its nucleotide variants keep their own
    fits, so the graphs keep their curves, and their class is `nonproductive`.
  - Where the amino-acid mutation count and parent link resolve, `pl7.app/vdj/kdLog10ParentDelta`
    and `pl7.app/vdj/kdLog10ParentDeltaSe` give each variant's log10 Kd,app against its parent's
    (labels "Kd,app vs Parent (log10)" and "Kd,app vs Parent SE (log10)"). Where the nucleotide
    mutation count and parent link (`[parentId, ntKey]`) resolve, the same two columns sit on the
    nucleotide axis: each nucleotide variant against the single unmutated nucleotide variant under
    its parent. A Failed variant has none; a variant censored at the edge of the titration reads
    from that edge. Parents that cannot be identified, and variants with several parents, are
    counted and named in the block's warnings.
  - Every output column on a variant axis carries a level step in its trace, as
    sort-seq-analysis does (type `milaboratories.synthetic-repertoire-profiler.level`, importance
    80): "NT Variants" on the nucleotide axis, "AA Variants" on the amino-acid axis. Outside this
    mode every column keeps the trace it had.
  - The Main table shows amino-acid rows in this mode. A new page, NT Variants, shows one row per
    nucleotide variant: its fit, its class and parent delta, the amino-acid variant it encodes (AA
    Variant Id, through the nt -> aa linker), and the result pool's columns on the nucleotide axis.
    Tabs show all nucleotide variants or only the synonymous variants of amino-acid variants with
    two or more fitted ones (`ntSynonymousOnly`, optional block data). The table state is new block
    data, added by a v7 -> v8 migration.
  - The graphs read each nucleotide variant's own fit, as in 2.5.0.

  Sample replicates:

  - The spread across replicates stays sample SD / sqrt(n), the standard error of a mean of
    independent replicates. The SE in `pl7.app/vdj/kdLog10Se` is the larger of that spread and the
    replicates' fit error, RMS(se_i) / sqrt(n), the rule synonymous variants follow, so replicates
    that agree more closely than their own fits allow do not rate Good. 2.5.0 reported the spread
    alone. `pl7.app/vdj/kdLog10Sd` and the replicate counts sit beside it, on every grain.
  - With synonymous variants as well, they combine inside each replicate first, and the amino-acid
    values then combine across replicates (Starr's order). The class, the stop rule and both parent
    deltas are decided once, on the final values. The amino-acid class and `kdOutOfRange`
    descriptions then state the rule across replicates.
  - A value combined from replicates, or from synonymous variants, rates Good only when more than
    half of them were measured inside the range. With both, the synonymous variants are counted
    inside each replicate: more than half of the replicates must each have more than half of their
    synonymous variants measured. One at an edge of the titration enters the mean at
    that edge, so two of them agree exactly and shrink the SE without measuring anything. A value
    with half or more at an edge rates Partial; its Kd,app is the one 2.5.0 reported. Under the mean
    bin fit, which reports no fit SE, the spread is the only evidence, so more than half of them
    must also be rated Good on their own; this holds across replicates too, where each amino-acid
    value carries the spread of its synonymous variants as its SE. With both, the rule
    counts the synonymous variants over the replicates that enter the mean. A bound is not a mean
    and keeps its class, with or without replicates.
  - A new page, Replicates, exists where the run had replicates. By Replicate (AA by Replicate in
    synonymous mode) lists the value each replicate contributed to the mean; in synonymous mode each
    row also carries the spread of the measured synonymous variants inside that replicate and their
    counts, and no R² or fit failure reason, which a combined value does not have. NT by Replicate
    (synonymous mode) lists every nucleotide variant's own fit in every
    replicate, beside its amino-acid variant (`per_replicate_nt.tsv`, `--out-per-replicate-nt`).
    Double-clicking a Main row still opens that variant's replicates.
  - The replicate settings that had no control in the UI are removed: the error model choice
    (`sd_over_sqrt_n`, `phillips`), the fit-error cutoff, the R² floor and the censoring majority.
    The `phillips` model reproduced that study's published SEM, which takes the standard error
    over the replicates and their own mean together: the population SD / sqrt(n + 1), 0.71 of the
    standard error of the mean for three replicates.
    The block runs one rule for each: sample SD / sqrt(n) floored by the fit error, no pre-combine
    filter, and a simple majority for an edge. A v7 -> v8 migration drops the stored values and
    the state of a per-replicate chart that no page shows (`graphStateReplicates`), the args carry the replicate column alone, the software loses the four `--replicate-*` option
    flags, and the `replicateFilteredCount` column is gone. The kind still accepts the four init
    params, so 2.5.x templates and scripts keep creating the block, and ignores them. A project
    with a replicate column goes stale once, because its args lose the two pinned values.
    Kd,app's description says that a censored replicate enters the mean at its edge. The
    per-replicate output drops its unread `kdState` column, and its -log10 Kd,app leaves an empty
    value out instead of storing it as NA.

  Fits and inputs:

  - Log fluorescence fit: a point whose reads all sit in an unsorted sample (no FACS bin) leaves
    before the points are counted. In 2.5.0 it failed the whole curve as `convergence_failure`. A fit
    rated Failed (by the R² or SE gate) writes no fitted curve, as under the mean bin fit. The reads
    of an unsorted sample count toward neither the read floor nor the fit weight of a point: only
    the sorted reads carry its signal. The mean bin fit is unchanged, as its output is pinned to
    2.4.0.
  - When the 0 M samples are unsorted, or none passes the read floor, each variant fits its own
    baseline, and the block says why in a warning.
  - A value censored at an edge of the range is rated by its wall (Below range, Above range) only
    under the log fluorescence fit. Under the mean bin fit it keeps the class its fits share, or
    Partial when they differ.
  - The (bin, concentration) grid check does not count an unsorted sample as a gate.
  - The FACS bin column can hold labels as well as integers. Picking a bin column lists its values,
    lowest fluorescence first, and the user drags them into order, as sort-seq-analysis does for its
    gates. Position is the gate's rank; no gate can be removed. Positive integers left in ascending
    order rank by their own value, so a gate number with no samples keeps its gap. A label column
    cannot run until its order is set. A bin value missing from the order fails the run and names
    the value; the settings name such values before the run and offer to add them. Two bins cannot
    share a rank. The order is new block data (`binOrder`, also a kind init param) and reaches the
    software as `--bin-ranks`. A block that never set an order reads the column's integers, as in
    2.5.0.
  - A log fluorescence fit with a FACS bin column but no gate fluorescence column shows a warning in
    the settings drawer: the fit then assumes that the gates are 0.65 decades apart, which can
    offset Kd,app.
  - `pl7.app/vdj/fitReads` ("Reads in fit", optional in the table) on both grains: the reads the
    fits used, every replicate included, from sorted samples at the non-zero concentrations that
    pass the read floor. The 0 M points set the baseline and are not counted, and a variant that was
    not fitted counts 0.
  - Each fitted curve is written at every concentration of the run, not only at the points its own
    fit saw, so a variant missing a point still has its curve there, and the mean curve across
    replicates reads every replicate at every concentration. Kd,app and the fits are those of
    2.5.0.
  - The synonymous-mode inputs are looked up only for a dataset whose variant axis may be
    nucleotide. Picking any other dataset records that (`synonymousInputs: false` in the args), and
    the block then does not wait on unrelated columns of the project.
  - A run whose abundance holds no read count for any selected sample stops with "the abundance
    column has no reads for any selected sample", instead of a polars type error.

  - An unsorted sample needs no sort fraction. The sort-fraction check reads sorted samples only;
    2.5.0 stopped the run when an unsorted sample had none.

  Pages and export:

  - The pages are Main, then the 2.5.0 graph pages unchanged: Titration Curves, Kd Distribution and
    Affinity vs Fit Quality. The NT Variants and Replicates tables have no page link; their values
    leave the block through Export raw data.
  - Export raw data on the Main page saves a ZIP of the tables the run wrote: the fits per variant,
    per amino-acid variant, per replicate and per nucleotide variant and replicate, the nt -> aa
    variant link, the reads each fit used, the observed signal and the fitted curves. Rows are keyed
    by the variant keys.
  - The Main table shows the upstream Mutations and mutation count columns and the count of
    synonymous nucleotide variants by default, and the variant class on request. Page titles and help text name the input: clonotype, variant or
    peptide. The tables list the Kd,app columns first, then the variant class, parent delta and
    counts, then the columns of other blocks.
  - The export carries R², the fit failure reason and `kdOutOfRange`, and under the mean bin
    fit the Hill coefficient, beside the affinity class, Kd,app, -log10 Kd,app and its SE, plus the
    parent delta and its SE where they resolve. 2.5.0 exported only the first four. The export
    also carries the spread: with replicates `kdLog10Sd` on every grain, without them
    `kdLog10SynonymSdAll` on the amino-acid grain, and `kdLog10SynonymSd` on the amino-acid grain
    either way. It carries the support behind a value too: `fitReads`, the replicate counts with
    replicates, and on the amino-acid grain the synonymous-variant counts. In synonymous mode it
    carries `variantClass` on both grains. The per-replicate values stay in the block's tables.
    The export does not carry the per-concentration signal (mean bin, fitted mean bin,
    concentration), which 2.4.0 exported.
  - The block builds on SDK model and ui-vue 1.88 and workflow-tengo 6.22.

### Patch Changes

- Updated dependencies [68af0de]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.6.0

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
