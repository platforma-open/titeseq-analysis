---
"@platforma-open/platforma-open.titeseq-analysis": patch
"@platforma-open/platforma-open.titeseq-analysis.model": patch
"@platforma-open/platforma-open.titeseq-analysis.workflow": patch
"@platforma-open/platforma-open.titeseq-analysis.ui": patch
"@platforma-open/platforma-open.titeseq-analysis.software": patch
---

UX review fixes:

- The block overview describes both fit equations, the five affinity classes, replicate
  combination and the -log10 Kd,app column.
- "Fail in-range fits of poor quality" and "Require reads at every concentration" have tooltips.
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
  replicates, and a hint above the table says so. The Open button no longer goes missing when the
  replicate outputs arrive after the table.
- With the log fluorescence fit, "FACS bin" is marked required and says why Run is disabled while it
  is empty. "FACS sort fraction" and "Gate fluorescence" show before a bin is picked.
- Log fluorescence fit: a Kd,app pinned at the edge of the range reports no R², as it reports no
  SE. With replicates, "Curve fit R²" is the mean over the in-range replicates, the ones Minimum R²
  checks. With the R² check on, no Good or Partial variant shows an R² below its minimum, and no
  bound shows an R² far below 0.
- "Equation" is now "Fit method equation", above the Fit method section.
