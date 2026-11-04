---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
---

Report one Kd,app per amino-acid variant on a nucleotide variant dataset, using synonymous nucleotide
variants as repeats.

- Turns on by itself when the dataset's variant axis is nucleotide and its nt -> aa linker
  (`pl7.app/repertoire/aaToNt`) resolves. Every other input runs as before.
- Each nucleotide variant is fitted on its own, exactly as before. Reads of two variants are never
  pooled. The fits are then combined per amino-acid variant: the mean of their log10 Kd,app. With
  two or more variants, the SE is the larger of their spread, SD / sqrt(n), and their own fit SEs
  combined, sqrt(sum of se_i^2) / n, where the equation reports a fit SE. With one variant, its own
  fit SE where the equation reports one. A variant censored at the edge of the titration enters at
  that edge; a Failed one is left out and counted.
- `pl7.app/vdj/kdLog10Sd` on the amino-acid axis: the SD of log10 Kd,app across the synonymous
  nucleotide variants in the mean, their spread alone with no fit-error floor (read-weighted under
  `reads`). Empty with fewer than two variants, for a bound and for a stop or failed variant.
- New setting `synonymousWeighting` (block data and the kind's init params; the drawer shows it only
  for a dataset that turns this mode on). Unset or `none`: the plain mean above. `reads`: each
  nucleotide variant is weighted by its total reads over the titration, in the mean, the spread and
  the fit-error term. Unset, the args and the command line are those of a block without the setting.
- Amino-acid columns sit on the linker's amino-acid axis and reuse the nucleotide names:
  `pl7.app/vdj/kd`, `pl7.app/vdj/kdNegLog10`, `pl7.app/vdj/kdLog10Se`, `pl7.app/vdj/affinityClass`,
  `pl7.app/vdj/kdOutOfRange`. The affinity class takes the equation's declared values and cutoff,
  as on the nucleotide grain. The export carries the same answer columns on both grains, plus the
  parent delta and its error where they resolve.
  `pl7.app/vdj/kdNtVariants` counts the nucleotide variants behind each value, split by
  `pl7.app/vdj/kdNtVariantState` (all, measured, censored low, censored high, failed).
- `pl7.app/vdj/variantClass` (String, discrete filter) on both grains, in Starr's strings: `wildtype`,
  `synonymous`, `stop`, `1 nonsynonymous`, `>1 nonsynonymous`. The amino-acid grain needs the
  amino-acid `pl7.app/repertoire/mutationCount`; the nucleotide grain needs that and the nucleotide
  one. `stop` needs the amino-acid `pl7.app/vdj/sequence/productive`; without it there is no stop
  class. Each input is used only when exactly one column matches it.
- A stop amino-acid variant has no Kd,app, SE or parent delta and is rated Failed; the run summary,
  the Fit log and the block's warnings count them. Its nucleotide variants keep their own fits, so
  the graphs keep their curves, and their class is `stop`. Starr nulls stop variants per barcode;
  this block's nucleotide grain keeps each raw fit.
- Where the amino-acid mutation count and parent link resolve, `pl7.app/vdj/kdLog10ParentDelta` and
  `pl7.app/vdj/kdLog10ParentDeltaSe` give each variant's log10 Kd,app against its parent's.
- Every output column on a variant axis carries a level step in its trace, as sort-seq-analysis does
  (type `milaboratories.synthetic-repertoire-profiler.level`, importance 80): "NT Variants" on the
  nucleotide axis, "AA Variants" on the amino-acid axis. A column on neither axis, and every column
  outside this mode, keeps the trace it had.
- Where the nucleotide mutation count and parent link (`[parentId, ntKey]`) resolve, the same two
  delta columns sit on the nucleotide axis: each nucleotide variant's own log10 Kd,app against the
  single unmutated nucleotide variant under its parent. A Failed variant has none; a variant
  censored at the edge of the titration reads from that edge. Parents that cannot be identified,
  and variants with several parents, are counted and named in the block's warnings, as at the
  amino-acid grain.
- The table shows amino-acid rows in this mode. Each nucleotide variant's own fit stays in the
  graphs, and a second page, NT Variants, shows one row per nucleotide variant: its fit, its class
  and parent delta, and the result pool's columns on the nucleotide axis. The page exists in this
  mode only. Its table state is new block data, added by a v6 -> v7 migration.
- Page titles and help text name the input: clonotype, variant or peptide.
- A run whose abundance holds no read count for any selected sample stops with "the abundance column
  has no reads for any selected sample", instead of a polars type error during normalization. This
  applies to every input kind.

With a sample replicate column as well, the copies combine inside each replicate first, and the
amino-acid values then combine across replicates (Starr's order). The class, the stop rule and both
parent deltas are decided once, on the final values. The amino-acid grain then carries the
replicates' SE, SD and counts, and the per-replicate drill-down shows amino-acid values.

Review fixes to the 2.5.0 fits and replicates:

- Log fluorescence fit: a point whose reads all sit in an unsorted sample (no FACS bin) leaves before
  the points are counted. It used to fail the whole curve as `convergence_failure`.
- When the 0 M samples are unsorted, or none passes the read floor, each variant fits its own
  baseline. The block now says why, as a warning, instead of logging "no 0 M samples". The software
  writes the run's warnings to `run_warnings.json`.
- A value censored at an edge of the range is rated by its wall (Below range, Above range) only
  under the log fluorescence fit. Under the mean bin fit it keeps the class its fits share, or
  Partial when they differ, so it never takes a value the column does not declare.
- The error across replicates is now population SD / sqrt(n + 1) (`pop_sd_over_sqrt_n_plus_1`)
  instead of sample SD / sqrt(n). On the block's own fits it lands closest to the published errors
  it was scored against. The setting value `phillips` is renamed to `pop_sd_over_sqrt_n_plus_1`. It is smaller, so more variants rate Good. A project with a replicate column
  goes stale on upgrade and gets the new error on its next run. The setting stays hidden.
- Every hidden Sample replicates setting is pinned to its default, not the error model alone, and
  none of them can disable Run.
- The (bin, concentration) grid check no longer counts an unsorted sample as a gate.
- Kd,app's description says that with replicates a censored replicate enters the mean at its edge.
- The per-replicate output drops its unread `kdState` column.

- The export again carries R², the fit failure reason and `kdOutOfRange`, and under the mean bin
  fit the Hill coefficient, beside the affinity class, Kd,app, -log10 Kd,app and its SE. 2.5.0
  exported only those four.
- A log fluorescence fit with a FACS bin column but no gate fluorescence column shows a warning in
  the settings drawer, which the user can close: the fit then assumes that the gates are 0.65
  decades apart, which can offset Kd,app. It is a new model output, `settingsWarnings`, and is not
  part of `validationWarnings`, so the results pages do not show it.

The export does not carry the per-concentration signal (mean bin, fitted mean bin, concentration),
which 2.4.0 exported. A downstream block that used one of those loses it on the next run of this
block.
