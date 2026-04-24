---
"@platforma-open/platforma-open.titeseq-analysis.workflow": patch
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
