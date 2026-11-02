---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
---

Choose the equation the fit uses, and import or export the fit settings as a file.

Two equations ship. `log_fluorescence` is new and is what a new block fits: a Hill fit on the log
gate fluorescence with the exponent fixed at 1, weighted least squares, tuned by a single plateau
floor. `hill_gate_index` is the fit the block shipped through 2.4.0, unchanged. Three presets name
the published analyses each reproduces — Adams 2016, Phillips 2021 (CR6261), and Starr 2020/2022 —
and write the settings that reproduce them.

A project created before this release is pinned to `hill_gate_index` by a migration that writes the
equation into the project's own data, so nothing about an existing fit moves. It also projects the
same workflow arguments it always projected, so the upgrade does not mark it stale: on the wire an
absent `equation` means `hill_gate_index`, which is deliberately not the default a new block takes.

The settings — and only the settings, never a column reference from another project — write out as
`titeseq-settings.json` and read back in. A picked file is read and checked in the model, which
reports the precise reason a file was refused; the settings are applied only when the scientist
presses Apply.

Separately, a clonotype whose signal never varies is now reported Failed rather than Partial. Such a
signal leaves the weighted R² undefined, and an undefined R² compares false against both quality
thresholds, so the row used to fall through to Partial carrying whatever the fitter returned — on a
flat bin mean that is a Kd,app of 28 M against a titration topping out at 1 µM. The row now reports
Failed with no Kd,app, exactly as a fit that never converged does. This is the one change to the
previous fit's behaviour: it moves only rows whose numbers were an artefact of the fit's amplitude
bound, and every clonotype in the regression corpus is unchanged.
