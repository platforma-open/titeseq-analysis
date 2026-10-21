---
'@platforma-open/platforma-open.titeseq-analysis.model': minor
'@platforma-open/platforma-open.titeseq-analysis.ui': minor
'@platforma-open/platforma-open.titeseq-analysis.workflow': minor
'@platforma-open/platforma-open.titeseq-analysis.kind': minor
'@platforma-open/platforma-open.titeseq-analysis': minor
---

Bin mode fits Kd,app with the Adams et al. 2016 model by default

New settings:
- **Sorted cell counts:** per-sample cells sorted into each bin. Required by the Adams model unless sort fractions are given or "fit without cell counts" is ticked.
- **Curve model section**, collapsed by default:
  - the model choice (Adams 2016 or linear Hill with free n);
  - bin width in decades;
  - minimum plateau;
  - the Good-class Kd,app uncertainty;
  - an optional 0 M baseline fluorescence.
- R², Hill-coefficient and hook settings now show only for the linear Hill model.
- The settings flag a missing 0 M control.

The summary gains two columns:
- **Kd,app uncertainty:** the posterior SD in decades.
- **Kd,app bound:** `below` or `above` when Kd,app is a limit of the titration.

Projects created before this version keep the linear Hill model they ran with. New blocks start on the Adams model.
