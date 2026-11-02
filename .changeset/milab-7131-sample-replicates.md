---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
---

Handle sample replicates in one block. Until now a block fitted one titration: samples sharing a
concentration were pooled into one point, and with the sort-fraction correction on the run stopped
with "sort_fraction values must sum to 1.0 per concentration", because two replicates sum to two.
Replicates had to be split into separate blocks.

- New optional setting `replicateColumnRef` (block data and the kind's init params): a per-sample
  String or integer column naming each sample's replicate. Unset, every output, every spec and
  every byte of the args and the command line is the one the block produced before. A column that
  cannot be used — missing, empty, or two values on one sample — warns and runs as one replicate.
- Each replicate is a titration of its own. The gate depth, the library depth and the sort-fraction
  sum are taken within (replicate, concentration), so the replicates no longer divide each other's
  reads, and the sort-fraction check passes on a dataset of several days' sorts.
- Each replicate is fitted alone; reads are never pooled across replicates. On a nucleotide variant
  dataset the order is Starr's: the synonymous copies combine inside each replicate first, then the
  replicates combine. The class, the stop rule and both parent deltas are decided once, on the
  combined value.
- The combined value is the plain mean of log10 Kd,app. A replicate whose Kd,app ran past a wall of
  the titration enters pinned at that wall, as Adams's and Phillips's code both do.
- New setting `replicateErrorModel`. `sd_over_sqrt_n` (the default) is the sample SD over the square
  root of the replicate count, falling back to the single fit's own standard error for one
  replicate. `phillips` is the population SD over the square root of the count plus one, which is 0
  for one replicate; it reproduces the standard error published with Phillips 2021.
  `pl7.app/vdj/kdLog10Sd` reports the sample SD across the replicates beside it.
- Two pre-combine filters, each a setting and each off unless set: `replicateFitErrorCutoff`
  (decades; Phillips 2021 published 1.0) and `replicateR2Floor` (Phillips 2021 published 0.8). Both
  apply to a replicate the fit measured inside the titrated range; a censored one is exempt, since
  the wall and not the fit set its value. A filtered replicate is left out of the mean and counted,
  and a variant whose every replicate is filtered is Failed with reason `all_replicates_filtered`.
- `replicateCensoredMajority` (default 0.5): when no replicate measured a Kd,app, one censoring side
  must hold at least this share of the replicates and exceed the other side for that wall to be
  reported; otherwise the variant is Failed with reason `conflicting_censoring`. This is the block's
  own rule — no published method defines one, and opposite walls cannot be averaged.
- `pl7.app/vdj/replicateCount` and its family — `replicateMeasuredCount`,
  `replicateCensoredBelowCount`, `replicateCensoredAboveCount`, `replicateFailedCount`,
  `replicateFilteredCount` — count what each variant's replicates produced, on both grains.
- `pl7.app/vdj/kdLog10Se` is now emitted under both equations when replicates are on: with several
  replicates the error is their spread, which does not depend on the equation.
- Two new frames on a (variant, replicate) axis, `replicatePf` and, in synonymous mode,
  `variantReplicatePf`, carry each replicate's own `pl7.app/vdj/kd`, `kdLog10Se`, `curveFitR2`,
  `affinityClass` and `kdReplicateState` (measured, censoredBelow, censoredAbove, failed,
  filtered). They are their own outputs: the extra axis would fan out every row of the main table.
  They stay out of the exported bundle.
- A new page, Replicates, shows those values as a table and as a chart that puts the replicates of
  one variant side by side. The page exists only when a run combined replicates. Its table and
  chart state is new block data, added by a v3 -> v4 migration that moves no argument.
- The signal frames gain the replicate on their key when a replicate column is set, so the
  Titration Curves graph draws each replicate's own curve instead of one pooled curve that belongs
  to no fit.
