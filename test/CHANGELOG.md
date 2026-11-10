# @platforma-open/platforma-open.titeseq-analysis.test

## 2.0.10

### Patch Changes

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

- Updated dependencies [68af0de]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.6.0

## 2.0.9

### Patch Changes

- @platforma-open/platforma-open.titeseq-analysis.model@2.5.2

## 2.0.8

### Patch Changes

- Updated dependencies [3e84ae6]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.5.1

## 2.0.7

### Patch Changes

- Updated dependencies [43eecfb]
- Updated dependencies [c6cb4c5]
- Updated dependencies [db4e101]
- Updated dependencies [db4e101]
- Updated dependencies [3ea39c8]
- Updated dependencies [7d02f2b]
- Updated dependencies [e776c4a]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.5.0

## 2.0.6

### Patch Changes

- Updated dependencies [03ca874]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.4.0

## 2.0.5

### Patch Changes

- e4569e0: Release every package through npm trusted publishing (OIDC) so each package's trusted publisher is used before it expires.
- Updated dependencies [e4569e0]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.2.1

## 2.0.4

### Patch Changes

- Updated dependencies [bccccbf]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.2.0

## 2.0.3

### Patch Changes

- Updated dependencies [d2ecf8c]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.1.2

## 2.0.2

### Patch Changes

- Updated dependencies [f42932b]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.1.1

## 2.0.1

### Patch Changes

- Updated dependencies [a3ad300]
  - @platforma-open/platforma-open.titeseq-analysis.model@2.1.0

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
