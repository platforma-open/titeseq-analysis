---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
---

Log fluorescence fit:

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
