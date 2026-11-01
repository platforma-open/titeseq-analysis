---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
---

Choose the equation the fit uses, with presets that reproduce published studies.

Two equations ship. `log_fluorescence` ("Log fluorescence fit") is new and is what a new block
fits: a Hill fit on the log gate fluorescence with the exponent fixed at 1, weighted least squares,
tuned by a single plateau floor. It reads an optional per-gate measured fluorescence column; without
it, the gates are assumed equally spaced. `hill_gate_index` ("Mean bin fit") is the fit the block
shipped through 2.4.0, unchanged. Three presets name the published analyses each reproduces —
Adams 2016, Phillips 2021 (CR6261), and Starr 2020/2022 — and write the settings that reproduce
them. Under `log_fluorescence` the Hill coefficient is a constant 1, so downstream blocks do not
receive it; the block's own table still shows it. Instead it outputs `pl7.app/vdj/kdLog10Se`, the
standard error of log10 Kd,app that rates its fits, and the Affinity vs Fit Quality page plots
Kd,app against it.

A project created before this release is pinned to `hill_gate_index` by a migration that writes the
equation into the project's own data, so nothing about an existing fit moves. It also projects the
same workflow arguments it always projected, so the upgrade does not mark it stale: on the wire an
absent `equation` means `hill_gate_index`, which is deliberately not the default a new block takes.

The fit now starts one worker per CPU the job is granted. It used to start one per host core, which
on a 96-core node ran 96 workers on a 4-CPU job. Sort fractions are now checked after the antigen
filter, so a multi-antigen dataset no longer fails the sum-to-1 check.

Separately, a clonotype whose signal never varies is now reported Failed rather than Partial. Such a
signal leaves the weighted R² undefined, and an undefined R² compares false against both quality
thresholds, so the row used to fall through to Partial carrying whatever the fitter returned — on a
flat bin mean that is a Kd,app of 28 M against a titration topping out at 1 µM. The row now reports
Failed with no Kd,app, exactly as a fit that never converged does. This is the one change to the
previous fit's behaviour: it moves only rows whose numbers were an artefact of the fit's amplitude
bound, and every clonotype in the regression corpus is unchanged.
