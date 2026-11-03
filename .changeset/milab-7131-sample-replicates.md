---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
---

Handle sample replicates in one block. Without a replicate column a block fits one titration:
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
