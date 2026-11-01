---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
---

Choose the fit equation, with presets that reproduce published studies.

- Two equations:
  - "Log fluorescence fit" (`log_fluorescence`), the new default: a Hill fit on log gate
    fluorescence with the exponent fixed at 1. Its one setting is the plateau floor. An optional
    column gives each gate's measured fluorescence. Without it, the gates are equally spaced.
    With 0 M samples, every clonotype shares one baseline, the read-weighted mean 0 M signal.
  - "Mean bin fit" (`hill_gate_index`): the 2.4.0 fit, unchanged.
- Presets set the values that Adams 2016, Phillips 2021 (CR6261) and Starr 2020/2022 used.
- Under `log_fluorescence` the block outputs `pl7.app/vdj/kdLog10Se`, the SE of log10 Kd,app. The
  Affinity vs Fit Quality page plots Kd,app against it, on a log y-axis. The Hill coefficient is a
  constant 1 and is not exported.
- A migration keeps existing projects on `hill_gate_index`. Their workflow arguments do not change,
  so the upgrade does not make them stale.
- The fit starts one worker per granted CPU, not one per host core.
- Sort fractions are checked after the antigen filter, so a multi-antigen dataset passes the
  sum-to-1 check.
- A clonotype with a constant signal is now Failed with no Kd,app. It was Partial, with a Kd,app set
  by the fit's amplitude bound.
