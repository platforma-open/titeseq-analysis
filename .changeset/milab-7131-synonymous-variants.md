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
- New setting `synonymousWeighting` (block data and the kind's init params; the drawer shows it only
  for a dataset that turns this mode on). Unset or `none`: the plain mean above. `reads`: each
  nucleotide variant is weighted by its total reads over the titration, in the mean, the spread and
  the fit-error term. Unset, the args and the command line are those of a block without the setting.
- Amino-acid columns sit on the linker's amino-acid axis and reuse the nucleotide names:
  `pl7.app/vdj/kd`, `pl7.app/vdj/kdLog10Se`, `pl7.app/vdj/affinityClass`, `pl7.app/vdj/kdOutOfRange`.
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
  export and in the graphs, and a second page, NT Variants, shows one row per nucleotide variant:
  its fit, its class and parent delta, and the result pool's columns on the nucleotide axis. The
  page exists in this mode only. Its table state is new block data, added by a v2 -> v3 migration.
- Page titles and help text name the input: clonotype, variant or peptide.
- A run whose abundance holds no read count for any selected sample stops with "the abundance column
  has no reads for any selected sample", instead of a polars type error during normalization. This
  applies to every input kind.
