# Concentration axis — spec deviation audit

**Date:** 2026-04-27
**Status:** Closed in this branch (`paulnewling/titeseq-user-feedback-fixes`).

## What the spec asks for

The affinity-profiling spec calls for a Float concentration axis on the
`meanBin` and `fittedMeanBin` PColumns:

- `docs/text/work/projects/affinity-profiling/pcolumn-spec.md:3`:
  *"mean bin and fitted mean bin have axes `[clonotypeKey][concentration:Float]`"*
- `docs/text/work/projects/affinity-profiling/pcolumn-spec.md:119`:
  *"The concentration axis is **Float** type to support log-scale rendering in
  Graph Maker."*
- `docs/text/work/projects/affinity-profiling/pcolumn-spec.md:123, 139`:
  schema blocks list `axes: [clonotypeKey][concentration:Float]` for each
  PColumn.
- `docs/text/work/projects/affinity-profiling/README.md:116` (R14):
  *"both with axes `[clonotypeKey][concentration:Float]`"*.

The spec also pins R14's parse-once canonical-string discipline for cross-layer
join stability:

> *"concentrations must be passed as their exact original string representation
> through the pipeline and converted to float only when used in arithmetic. The
> concentration axis key written to the PColumn must use this canonical string,
> not a value that was parsed and re-serialized."*

## Why the spec is unimplementable in the current SDK

The workflow-tengo `pt` layer (which `xsv.importFile` ultimately delegates to)
gates axis types to `Int | Long | String`:

- `core/platforma/sdk/workflow-tengo/src/pt/util.lib.tengo:352` —
  `` `type`: `string,regex=Int|Long|String` `` on the axis schema.
- `core/platforma/sdk/workflow-tengo/src/pt/import-dir.tpl.tengo:33-38` —
  side-by-side asymmetry showing the constraint is deliberate, not an oversight:

  ```
  `axes`: [{
      `id`: `string`,
      `type`: `string,regex=Int|Long|String`            ← axis types
  }],
  `column`: {
      `id`: `string`,
      `type`: `string,regex=Int|Long|Float|Double|String` ← column types
  }
  ```

  `Float` and `Double` are explicitly excluded from axis types and explicitly
  included in column types. Two adjacent regex declarations, in the same file,
  make the asymmetry intentional.

Empirically verified: a Float-axis attempt during this branch's iteration
produced the runtime error
`Invalid params structure: axes[1].spec.type: value "Double" does not conform
regex "Int|Long|String"` (trace through `pframes.xsv.importFile` →
`pframes.xsv-import-pt` → `pt:1329` → `pt.util:37`).

A search across the workspace confirmed **no block emits a Float or Double
axis** — every `Float` / `Double` use in `blocks/*/workflow/` is on a column,
never an axis. The spec is currently aspirational.

## Prior workaround (the molar Long aM axis)

The original implementation chose `axes [clonotypeKey, concentrationStr:String,
concentrationAM:Long]` where `concentrationAM = round(value * 1e18)` (Molar →
attomolar). Long was the only numeric axis type the SDK accepts; attomolar
kept typical TiteSeq concentrations (`10⁻¹³ … 10⁻⁶ M`) inside int64 range.

This baked in two problems:

1. **R2 violation.** `× 1e18` is a unit assumption (input is molar). R2 says
   the block treats values as dimensionless floats. A user entering `100` (nM)
   hit the `MAX_CONCENTRATION_M ≈ 9.2` ceiling with a misleading "molar vs
   molal" hint, even though `100` nM is well within physical range.
2. **Graph X-axis ticks at `10⁶ … 10¹²`** (aM magnitudes) are unreadable to a
   user thinking in molar / nM / µM. Required mental conversion via 10⁻¹⁸.

A defensive dual-import pattern (`*Internal` 3-axis variant + `*Export` 2-axis
variant with `pl7.app/hideDataFromGraphs`) was added to keep the Long axis out
of downstream blocks' graph pickers. See
`docs/investigations/cid-conflict-investigation.md` for the CID-conflict risk
this pattern introduced.

## What this branch ships instead

A **String axis** as the canonical concentration axis:

```
meanBin            axes: [clonotypeKey, concentration:String]   valueType: Double
fittedMeanBin      axes: [clonotypeKey, concentration:String]   valueType: Double
```

- The String axis IS the canonical join key. R14's parse-once discipline is
  honored at the byte level — the Tengo workflow wraps `concentrationStr`
  directly as a String axis on the output PColumns; no float parse, no
  re-serialization.
- Graph Maker binds X to the `concentration:String` axis directly. The
  Titration Curves scatter renders the X-axis categorically, ordered by
  parsed numeric value (Graph Maker's String-axis sort is numeric-aware).
- Single import per TSV. The same pcolumn resource is reused in `signalPf` and
  `exportPf` — no `*Internal`/`*Export` variants, no `hideDataFromGraphs` hack,
  closes the dual-import CID-conflict pattern.
- No `× 1e18` anywhere. No `MAX_CONCENTRATION_M` ceiling. The pipeline is
  unit-agnostic — any input unit (M, nM, µM, …) works.

### Tradeoff: categorical X-spacing vs. additionalCurves

An earlier iteration of this branch added a `concentrationValue` Double sidecar
PColumn so X could be bound to a numeric source and rendered on a true log
scale. That pattern broke the `additionalCurves` slot in Graph Maker's
scatter chart: `pf-plots@1.2.0`'s `checkSourceBySpec` gate requires both X and
Y to resolve to **axis specs**, and silently greys out additional curves when
either one is bound to a PColumn spec. With X bound to a PColumn, the Hill
curve overlay was unreachable.

The fix on the block side: bind X to the `concentration:String` axis directly
and drop the sidecar. We accept categorical (equal-width) X-tick spacing
instead of true log-scale spacing in exchange for the curve overlay. Tick
labels still render at the user's input concentration values (`5e-6`, `1e-6`,
…), so visually a user reads the same x positions — only the spacing between
them changes.

The proper fix lives in the SDK / pf-plots: see
`docs/investigations/sdk-pr-additional-curves-pcolumn-x.md` for the planned
upstream PR. When that ships, this block can re-introduce the sidecar
PColumn for true log spacing.

## Spec alignment scorecard

| Spec requirement | Prior implementation | This branch |
|--|--|--|
| R2: "values are dimensionless floats" | violated (× 1e18 = molar assumption) | honored — no unit assumption |
| R14: `axes [clonotypeKey][concentration:Float]` | unimplementable in SDK | closest available — `[clonotypeKey][concentration:String]` |
| R14: parse-once canonical string preservation | preserved (parallel `concentrationStr` axis) | preserved — String axis IS the canonical key |
| `pcolumn-spec.md:119` "Float type to support log-scale rendering" | Long aM (workaround) | partial — categorical X-axis ordered by numeric value, not true log spacing (blocked on `additionalCurves` SDK bug; see `sdk-pr-additional-curves-pcolumn-x.md`) |
| R15: "x = concentration (log scale)" | renders aM ticks (`10⁶…10¹²`) | renders user's input values directly, categorical spacing |
| R2: Kd's `pl7.app/unit` annotation from concentration column label | works | unchanged |
| Misleading `MAX_CONCENTRATION_M` rejection | present | gone |
| Hill curve overlay (`additionalCurves`) | works | works (depends on X binding to axis) |

## Forward path

Two upstream gaps keep this block off the literal spec:

1. **Axis-type regex.** A future SDK PR widening axis types to include
   `Float | Double` would make the spec literally implementable, at which
   point this block could put floats on the axis directly and drop the
   String-axis workaround.
2. **`additionalCurves` PColumn-X support.** The pf-plots
   `checkSourceBySpec` gate currently rejects PColumn specs as the X source
   for the additionalCurves overlay slot, even though Graph Maker handles
   PColumn-bound X correctly elsewhere. Once the gate accepts PColumn
   sources, this block can re-introduce a Double sidecar PColumn alongside
   the String axis to get true log-scale X spacing without losing the Hill
   curve overlay. See `docs/investigations/sdk-pr-additional-curves-pcolumn-x.md`.

Per the team's spec-doc workflow (memory `feedback_spec_docs.md`: do not edit
spec docs), the correction lives here — code comments at the deviation site
(`workflow/src/main.tpl.tengo`'s `concAxisSpec` block) plus this audit.

## References

- Workflow change: `workflow/src/main.tpl.tengo` (`concAxisSpec`, single-import
  pattern).
- Python pipeline: `software/src/output_build.py` (`build_mean_bin_frame` —
  outputs canonical `concentrationStr` only).
- UI: `ui/src/pages/TitrationCurvesPage.vue` (X bound to the
  `concentration:String` axis).
- SDK follow-up: `docs/investigations/sdk-pr-additional-curves-pcolumn-x.md`
  (proposes the pf-plots fix that would let X bind to a Double sidecar
  PColumn while keeping the curve overlay — would re-enable true log-scale
  X spacing here).
- Earlier related fix: `docs/investigations/cid-conflict-investigation.md`
  (params canonical-encoding + dual-import discussion).
