---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": patch
---

Statistics and review fixes, and removal of two pages that had no link.

Statistics and review fixes:

- Parent delta with sample replicates. The block takes the delta inside each replicate that
  measured both the variant and its parent inside the range. It then combines the deltas as it
  combines replicate values. A shift that a replicate gives every value cancels. A replicate that
  measured only one of the two is left out of the delta, so the delta can differ from the
  difference of the two Kd,app values beside it by a shift between replicates. Where a replicate
  that has both holds the variant or the parent at an edge of the range, or no replicate measured
  both, the delta compares the combined values and has no SE.
- New columns `pl7.app/vdj/kdParent` and `pl7.app/vdj/kdLog10ParentSe` hold the parent's own Kd,app
  and SE. They are optional table columns, and both grains export them.
- A run with sample replicates logs how many parent deltas compare the combined values, and how
  many of those lie between two values inside the range. The description of the delta's SE column
  says when it is empty.
- Log fluorescence fit: R² counts the variant's own 0 M points, from every 0 M sample. Each point
  counts at its read weight. Together they weigh at most the variant's heaviest fit point. A tight
  binder near its plateau passes the R² gate. A deeper 0 M sample, or one more 0 M sample, does not
  raise R² above that cap. A variant with no 0 M point past the read floor counts the shared
  baseline as its 0 M point, at the weight of its heaviest fit point. The R² gate removes a curve
  that does not rise from the baseline; the SE of log10 Kd,app rates how well a strong binder is
  fitted.
- Log fluorescence fit: a variant that is flat and bright at every point, its own 0 M point
  included, is rated Failed with the reason `bright_without_antigen`. It was rated Below range, as
  a saturated tight binder is. A saturated binder with a 0 M point rises from it and is fitted. A
  variant with no 0 M point of its own past the read floor is rated as before.
- Mean bin fit with a sort fraction: Starr 2022's split rule. A concentration has no mean bin and
  leaves the fit when two gates at least two apart each hold more than 40% of the variant's cells.
  Without a sort fraction the rule is off, because shares of reads over gate depth split a single
  peak when the gates are sequenced unevenly.
- Both fits: a variant fails as `insufficient_points`, not `insufficient_reads`, when it has
  reads past the read floor at a concentration but no point is left. A point is left out by the
  split rule, or when no gate that holds its reads carries weight.
- `argsIssues` disables Run and names every problem at once. The args no longer throw. Valid
  settings send the same args as before.
- The settings drawer shows the message for a sort fraction without a bin column. The desktop
  names no reason beside a disabled Run.
- Model outputs read each bin column once: the output `binColumns` replaces `binOptions` and
  `binFactsByRef`. They skip the result pool scan for a dataset whose
  variant axis is not nucleotide.
- A project keeps its results until it is run again with changed settings. The new rules apply
  from that run.

Pages removed:

- The NT Variants and Replicates pages had a route and no link. They are removed with their tables
  (`ntVariantTable`, `replicateTable`, `replicateNtTable`) and the two frames only they read
  (`ntVariantPf`, `replicateNtPf`).
- A variant's replicates still open from its row on Main.
- Every nucleotide variant's and replicate's values stay in the raw data download. The nucleotide
  grain's class and parent delta join it as `nt_variant_annotations.tsv`.
- The table states of the removed pages stay in the block data, unread. No migration is needed.
- `nt_to_aa_variant.tsv` in the raw data download pairs each nucleotide variant with its
  amino-acid variant. Its `aaVariant` label column fed the NT Variants page alone and is removed;
  both variants are now named as in every raw data file (below).
- With sample replicates, the synonymous combination inside each replicate no longer runs the
  class and parent delta annotation. The combined values decide both. Results are unchanged.

Affinity vs Fit Quality and downloads:

- The `Kd,app (plot position)` and `Hill coefficient (plot position)` columns are removed. They
  parked a missing Kd,app or Hill coefficient at a fixed value. Affinity vs Fit Quality plots
  Kd,app against the SE of log10 Kd,app (log fluorescence fit) or the Hill coefficient (mean bin
  fit), grouped by affinity class. Under the mean bin fit it shows every class, not only Failed.
  A fit with no Kd,app is not drawn; the table lists it with its failure reason.
- A saved Affinity vs Fit Quality chart that drew a plot position starts over from the default.
- The tables lose their Export button. The raw data download is the one export.
- Every file in the raw data download names each variant by its Variant Id and sequence, in place
  of its hash key. `nt_to_aa_variant.tsv` names both grains, each prefixed NT or AA. The sequence
  is the upstream amplicon sequence column in the variant's own alphabet; it is empty when the
  upstream has none. A variant with no label column keeps its key as its Variant Id.
- Every column in the raw data download is headed by the label the results table shows for it
  (`Kd,app`, `Affinity class`, ...), and the variant columns by the upstream columns' labels. A
  column the block does not import keeps its own name.
