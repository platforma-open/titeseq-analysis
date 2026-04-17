# Spec review — titeseq-analysis block implementation

Review of `model/src/index.ts`, `workflow/src/main.tpl.tengo`, and `ui/src/**` against `docs/text/work/projects/affinity-profiling/README.md` (R1–R23) and `pcolumn-spec.md`. Python software in `software/` is treated as complete. Build passes: `pnpm run build:dev` is green.

## R1–R5 Inputs

| Req | Status | Location | Notes |
|-----|--------|----------|-------|
| R1 abundance `[sampleId][clonotypeKey]` | ✅ | `model/src/index.ts:142-160`; `workflow/src/main.tpl.tengo:37-38,71-73` | Supports both `clonotypeKey` and `scClonotypeKey` |
| R2 concentration + unit propagation | ✅ | `model/src/index.ts:164-172,196-200`; `workflow:61-64,150-160` | Unit copied into `pl7.app/unit` on `kd` |
| R3 bin column (optional) | ✅ | `model/src/index.ts:173-183`; `workflow:41-42,58-80` | Integer `[sampleId]` filter |
| R4 antigen + targetAntigen | ✅ | `model/src/index.ts:184-193`; `workflow:44-46,117-122` | Model blocks run if antigenRef set without targetAntigen; workflow forwards both to Python |
| R5 1 sampleId = 1 value | ⚠️ | Deferred to Python | Validation happens in `io_layer.py` (out of scope per user note) |

## R6–R12 Curve Fitting & Classification

All delegated to Python (`software/src/...`); block surfaces results via `per_clonotype.tsv` and `mean_bin.tsv`/`fitted_mean_bin.tsv`. User said to treat Python as complete.

**Thresholds piped through `params.json`** (`workflow:86-97`): `min_reads_per_concentration`, `min_concentration_points`, `r2_threshold_good`, `r2_threshold_failed`, `n_min`, `n_max`, `hook_effect_threshold_bin`, `hook_effect_threshold_no_bin`, `hook_effect_min_reads`. Matches R10/R12 defaults in the spec.

## R13–R14b Outputs

| Req | Status | Location | Notes |
|-----|--------|----------|-------|
| R13 summary cols `[clonotypeKey]`; meanBin `[clonotypeKey][concentration]` | ✅ | `workflow:165-257, 262-309` | Correct axes for all 8 output columns |
| R14 meanBin + fittedMeanBin w/ Float conc axis; c=0 excluded | ✅ | `workflow:129-140, 262-309` | Float axis; c=0 filtering in Python |
| R14b kdOutOfRange flag | ✅ | `workflow:241-253` | `pl7.app/vdj/kdOutOfRange` Boolean |

## PColumn spec match (pcolumn-spec.md)

All 8 columns match exactly:

| Column | Name | ValueType | Score | Visibility | OrderPri | Status |
|--------|------|-----------|-------|------------|----------|--------|
| Affinity class | `pl7.app/vdj/affinityClass` | String | true | default | 202 | ✅ |
| K_D,app | `pl7.app/vdj/kd` | Double | true | default | 201 | ✅ (unit injected) |
| Curve fit R² | `pl7.app/vdj/curveFitR2` | Double | false | default | 199 | ✅ (renamed from Python `r2` via column mapping) |
| Hill coefficient | `pl7.app/vdj/hillCoefficient` | Double | false | default | 198 | ✅ |
| Fit failure reason | `pl7.app/vdj/fitFailureReason` | String | false | hidden | 197 | ✅ |
| K_D,app out of range | `pl7.app/vdj/kdOutOfRange` | Boolean | false | default | 196 | ✅ |
| Mean bin | `pl7.app/vdj/meanBin` | Double | false | hidden | — | ✅ |
| Fitted mean bin | `pl7.app/vdj/fittedMeanBin` | Double | false | hidden | — | ✅ |

Affinity class discrete metadata (`pl7.app/isDiscreteFilter`, `pl7.app/discreteValues`, `pl7.app/score/rankingOrder`, `pl7.app/score/defaultCutoff`) all present at `workflow:217-220`.

## R15–R19b Visualizations

| Req | Status | Location | Notes |
|-----|--------|----------|-------|
| R15 Titration Curves — scatter + curve, facet clonotypeKey, log-x | ✅ | `ui/src/pages/TitrationCurvesPage.vue` | Default options set x=concentration axis, y=meanBin, grouping=clonotypeKey. Both meanBin + fittedMeanBin live in `signalPf` (workflow:331-338), so Graph Maker can render both in one PFrame. Log scale via `graphStateTitrationCurves.axesSettings.axisX.scale = "log"` (model:92-95). **Note: curveOverlay slot (R22) is cross-project in visualization-api; rendering currently reads both columns but R22 enables true overlay.** |
| R16 K_D histogram on log x | ✅ | `ui/src/pages/KDDistributionPage.vue` | `chart-type="histogram"`, filter `affinityClass ∈ {Good, Partial}` (lines 24-29). Log x scale in UiState default (model:104-107) |
| R17 Affinity vs fit, all clonotypes incl. Failed | ✅ | `ui/src/pages/AffinityVsFitPage.vue` | x=kd, y=hillCoefficient, grouping=affinityClass. No filter = includes Failed. `fitFailureReason` exists as separate column for per-failure encoding. Log x in UiState (model:119-121) |
| R18 Table tab, sortable | ✅ | `ui/src/pages/TablePage.vue` | `PlAgDataTableV2` with `usePlDataTableSettingsV2({ model: () => app.model.outputs.summaryTable })` |
| R19 No workflow rerun on param change | ✅ | GraphMaker state + tableState | All viz state is UI-only; Graph Maker and PlAgDataTable recompute locally |
| R19b no-bin mode warning + relabel | ✅ | `TitrationCurvesPage.vue` PlAlert + `workflow:142-148` | Label swap "Mean bin" → "Clonotype frequency" at the spec level; persistent warning banner on TitrationCurvesPage |

## R20–R21 Lead Selection integration

| Req | Status | Location | Notes |
|-----|--------|----------|-------|
| R20 Discovery annotations | ✅ | `workflow:165-257` | All score columns carry `pl7.app/isScore: "true"` + `pl7.app/vdj/` prefix |
| R21 Exports visible to Lead Selection | ✅ | `workflow:341-351, 359-361` | `exports.pf` contains summary + signal columns merged |

## R22–R23 Graph Maker co-development

**Deferred** — these are visualization-api work, not block work. Block outputs the data in the shape the curve overlay slot expects (`fittedMeanBin` aligned to `meanBin` axes). ✅

## Block structure & plumbing

- **Model validation (R2, R4)**: `argsValid` blocks run until antigenRef without targetAntigen, threshold inversions, out-of-range params. Non-blocking warnings surfaced via `validationWarnings` output (spaces-in-unit warning, targetAntigen-without-ref warning).
- **Title + subtitle**: R2's unit label appears via `title()` when abundance ref is set (model:311-321).
- **Settings drawer auto-open**: `PageHeader.vue` opens settings if `abundanceRef === undefined`.
- **Log viewer**: `PageHeader.vue` wraps `PlLogView` with `app.model.outputs.logHandle` (from `fit.getStdoutStream()`).
- **isRunning / isEmpty**: both derived from workflow output status; `isEmpty` drives empty-state banner in TitrationCurvesPage.
- **Trace injection**: `pSpec.makeTrace` + `trace.inject` applied to all output specs (workflow:318-322, 325-351). `splitDataAndSpec: true` means only ~8 small spec resources are rewritten; data resources pass through.

## Known gaps & follow-ups

1. **R5 validation is Python-side only.** The model does no per-sample uniqueness check (would need a dedicated workflow step). Python's `io_layer.py` rejects conflicting sampleId rows — consistent with spec "Reject the block with a user-facing error" but gated by Python rather than Tengo.
2. **R15 curve overlay is a single-source scatter today.** Without the visualization-api `curveOverlay` slot (R22), Graph Maker renders both `meanBin` and `fittedMeanBin` as dots. The data shape is already correct for the overlay slot when it lands — no block-side changes needed.
3. **SDK version warnings.** `pnpm install` emits benign peer warnings for `vite`, `@milaboratories/pl-model-common`, `vitest`. Build is clean. Before opening a PR, run `npm view @platforma-sdk/block-tools version` and bump `pnpm-workspace.yaml` catalog if stale.
4. **Tag for local filtering.** Add `"paul-newling"` to `block/package.json` `block.meta.tags` as an unstaged local modification per project convention.
5. **Changeset.** Need `.changeset/<name>.md` before PR.

## Result

Implementation hits R1–R21 block-side requirements. Build: ✅. All pcolumn-spec.md annotations match byte-for-byte. R22/R23 are cross-project work outside this block. R5 validation is pushed into Python per user's "treat Python as complete" guidance.
