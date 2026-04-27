---
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis": minor
---

- Switch `params.json` encoding from the generic `json` module to
  `@platforma-sdk/workflow-tengo:canonical` so repeated runs produce
  byte-identical input. Eliminates a transient `CIDConflictError` path.
- Bump `@platforma-sdk/*` catalog entries (workflow-tengo 5.15.0,
  block-tools 2.7.13, model/ui-vue/test 1.67.0, tengo-builder 2.5.14).
- Expand `validationWarnings` to surface backend output errors
  (CID conflicts, exec failures) and mirror every `.args()` invariant as
  a severity:"error" issue — a disabled Run button now always has a
  visible reason.
- Render validation alerts on every page (not just the settings drawer)
  via the shared `TiteseqPage` shell.
- Enforce integer domains on `minReadsPerConcentration`,
  `minConcentrationPoints`, and `hookEffectMinReads` with inline
  "Must be a whole number" errors, and mark them `required` so the
  red-asterisk styling matches the ref-picker inputs.
- Enable CSV export on the Clonotype Fit Results table via the
  SDK-native `show-export-button` prop, and add a new "Mean Bin Data"
  section that renders `signalPf` as a `PlAgDataTableV2` with the same
  export button. Together these cover the per-clonotype Hill fit
  outputs and the per-concentration mean-bin/fitted-mean-bin data the
  curves are derived from.
- Realign the concentration axis to spec. Previously the workflow emitted
  a Long attomolar axis (`× 1e18`) that baked in a hidden molar
  assumption, violating R2 ("values are dimensionless floats") and
  rendering graph X-axis ticks at `10^6 … 10^12` (alien aM magnitudes).
  Now uses a single String axis (canonical concentration string is the
  join key, per R14). Graph Maker binds X to that axis directly, so
  X-axis ticks render at the user's actual input concentrations. Drops
  the misleading `MAX_CONCENTRATION_M ≈ 9.2` rejection — any input unit
  / magnitude now works. Single import per TSV (no more dual-pcolumn
  `_Internal`/`_Export` variants), which also closes the dual-import
  CID-conflict pattern as a side effect. Tradeoff: X-axis spacing is
  categorical (ticks ordered by numeric value, equal width between
  them), not true log-scale spacing. This is required so the Hill
  curve overlay (`additionalCurves`) un-greys; pf-plots'
  `checkSourceBySpec` rejects PColumn-bound X sources for the curve
  slot. Spec calls for `concentration:Float` but the SDK gates axis
  types to `Int|Long|String` — see
  `docs/investigations/concentration-axis-spec-realignment.md` for the
  deviation rationale and `docs/investigations/sdk-pr-additional-curves-pcolumn-x.md`
  for the upstream fix that would restore true log spacing.
