---
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

Fix wrong results and failed runs found in the 2026-10-20 bug hunt

- A clonotype with a flat signal is Failed (`low_r2`), not Partial with an invented Kd,app.
- The linear Hill fit's Kd,app range follows the data's concentration unit; with nM concentrations, weak binders no longer read exactly 1000.
- An unsorted 0 M control no longer pulls the mean bin at 0 M, and with it the baseline, to 0.
- A sort fraction combined with an antigen column no longer fails the run: sums are checked within the target antigen.
- A pre-sort library sample without a sort fraction or cell count no longer fails the run.
- With too few concentrations for the Hill fit's parameters, clonotypes are `insufficient_points`, not `convergence_failure`.
