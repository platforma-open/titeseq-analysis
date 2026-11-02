---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
---

Report one Kd,app per amino-acid variant on a nucleotide variant dataset, using synonymous nucleotide
variants as repeats.

- Turns on by itself when the dataset's variant axis is nucleotide and its nt -> aa linker
  (`pl7.app/repertoire/aaToNt`) resolves. There is no setting. Every other input runs as before.
- Each nucleotide variant is fitted on its own, exactly as before. Reads of two variants are never
  pooled. The fits are then combined per amino-acid variant: the mean of their log10 Kd,app, with
  SE = SD / sqrt(n). With one variant, its own fit SE where the equation reports one. A variant
  censored at the edge of the titration enters at that edge; a Failed one is left out and counted.
- Amino-acid columns sit on the linker's amino-acid axis and reuse the nucleotide names:
  `pl7.app/vdj/kd`, `pl7.app/vdj/kdLog10Se`, `pl7.app/vdj/affinityClass`, `pl7.app/vdj/kdOutOfRange`.
  `pl7.app/vdj/kdNtVariants` counts the nucleotide variants behind each value, split by
  `pl7.app/vdj/kdNtVariantState` (all, measured, censored low, censored high, failed).
- Where the amino-acid mutation count and parent link resolve, `pl7.app/vdj/kdLog10ParentDelta` and
  `pl7.app/vdj/kdLog10ParentDeltaSe` give each variant's log10 Kd,app against its parent's.
- The table shows amino-acid rows in this mode. Each nucleotide variant's own fit stays in the
  export and in the graphs.
- Page titles and help text name the input: clonotype, variant or peptide.
