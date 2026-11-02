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
- Where the amino-acid mutation count and parent link resolve, `pl7.app/vdj/kdLog10ParentDelta` and
  `pl7.app/vdj/kdLog10ParentDeltaSe` give each variant's log10 Kd,app against its parent's.
- Every output column on a variant axis carries a level step in its trace, as sort-seq-analysis does
  (type `milaboratories.synthetic-repertoire-profiler.level`, importance 80): "NT Variants" on the
  nucleotide axis, "AA Variants" on the amino-acid axis. A column on neither axis, and every column
  outside this mode, keeps the trace it had.
- The table shows amino-acid rows in this mode. Each nucleotide variant's own fit stays in the
  export and in the graphs.
- Page titles and help text name the input: clonotype, variant or peptide.
