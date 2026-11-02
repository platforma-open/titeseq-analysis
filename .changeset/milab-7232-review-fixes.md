---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
---

Log fluorescence fit:

- A flat, dim curve rates `Non-binder` instead of `Partial`. Lead Selection's default affinity
  cutoff keeps Good and Partial, so non-binders no longer pass through it as candidates. A curve
  that rose but did not reach half-maximum inside the range stays `Partial`.
- A new column, "Kd,app vs titration", says whether Kd,app sits below, inside or above the titrated
  range. "Kd,app out of range" carries the same fact without the side and moves to the optional
  columns.
- An optional R² gate: an in-range fit with R² below the threshold rates `Failed`, `low_r2`. Off
  by default. The Phillips 2021 preset turns it on at 0.8, as its authors filter. Curve fit R² is
  hidden while the gate is off, and is empty instead of NaN for a signal with no variance.
- Gate fluorescence above 7 log10 is refused as linear fluorescence, with an error that names the
  column. It used to move Kd,app by decades, or overflow the fit.
- With no FACS bin column, the reason Run is disabled names the fix: choose the Mean bin fit.

Both fits:

- The missing (bin, concentration) warning counts the clonotypes of the selected antigen only, as a
  share of them, and no longer prints internal clonotype keys.
- The Kd,app unit is the text in the concentration label's trailing parentheses: a column labelled
  "Antigen concentration (M)" gives "Kd,app (M)". Kd,app carries a description of what "apparent"
  means.

The mean bin fit's columns and classes are unchanged.
