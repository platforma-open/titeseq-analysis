# SDK fix plan — `additionalCurves` should accept PColumn X/Y sources

## Read this first — version skew

The bug being fixed here lives in **published `@milaboratories/pf-plots@1.2.0`**, which is what `titeseq-analysis` (and likely other downstream blocks) depends on. The local clone of `milaboratory/visualizations` is at version **`1.1.66`** on `main` (latest `CHANGELOG.md` entry there is `1.1.66`). **The source for `1.2.0` is NOT in the local repo** — `git log -S 'additionalCurves'` on `core/visualizations` returns zero matches across `--all`.

What this means for whoever picks this up:

- **Don't expect to find `additionalCurves` or `checkSourceBySpec` in `core/visualizations/packages/pf-plots/src/`** — the local clone genuinely does not contain them. A previous version of this plan claimed otherwise; that was wrong.
- The bug is real. It's verifiable by inspecting `node_modules/.pnpm/@milaboratories+pf-plots@1.2.0/.../dist/controllers/controllersByChartType/scatterplot.js` and its `.js.map`. The minified `scatterplot.js` contains `additionalCurves` 6 times; the source map's `sourcesContent` contains both `additionalCurves` and `checkSourceBySpec` in the original TypeScript form.
- Conclusion: `1.2.0` was published from a build/source that was never committed to the public `main` of `milaboratory/visualizations`. Either an unpushed branch, a private fork, or a manual local build.

## What to do first

Before writing any code, resolve the source question:

1. **Ask the pf-plots maintainer** (whoever publishes `@milaboratories/pf-plots`) where the `1.2.0` source lives. If there's a branch or fork to clone, get it and proceed normally.
2. **If the source is genuinely missing** (lost, unpushed, or otherwise unrecoverable), treat this as a **feature add** against the local `1.1.66`: add `additionalCurves` to scatterplot from scratch with the corrected semantics described below, and publish a new version (`1.2.1` or `1.3.0`) as a successor to the orphaned `1.2.0`.

The rest of this doc assumes the second scenario (feature add against `1.1.66`), since that's the worst case and produces a self-contained PR. If the actual `1.2.0` source surfaces, the same semantic fix applies — just patched into the existing controller instead of added from scratch.

## Why the change is needed (titeseq-analysis use case)

The `affinity-profiling` spec ([excerpts in `docs/investigations/concentration-axis-spec-realignment.md`]) calls for a Float concentration axis on `meanBin` and `fittedMeanBin`. The workflow-tengo SDK's `pt` layer rejects Float on axes (`core/platforma/sdk/workflow-tengo/src/pt/util.lib.tengo:352` — axis types regex `Int|Long|String`). The block's workaround:

```
meanBin            axes: [clonotypeKey, concentration:String]   valueType: Double
fittedMeanBin      axes: [clonotypeKey, concentration:String]   valueType: Double
concentrationValue axes:                       [concentration:String]   valueType: Double  ← sidecar
```

In Graph Maker, the user binds:
- **X** = `concentrationValue` (a PColumn — Double values let Graph Maker render true log scale)
- **Y** = `meanBin` (a PColumn)
- **additionalCurves** = `fittedMeanBin` (the Hill curve overlay)

All three share the `concentration:String` axis. The dot rendering works correctly (verified visually — X-axis ticks at `5×10⁻⁶`, `1×10⁻⁶`, etc. on log scale). The Hill curve overlay does not — the source picker for `additionalCurves` is universally greyed out in the desktop UI.

What I observed in the published `1.2.0` source map (the bug):

```ts
additionalCurves: new ComponentController<ScatterplotUIState, 'additionalCurves'>({
  componentName: 'additionalCurves',
  allowedTypes: ['Int', 'Long', 'Double', 'Float'],
  strictlyDependsOnParents: true,
  parentComponents: ['x', 'y'],
  dependsOn: ['x', 'y'],
  settings: { columnsAllowed: true, multipleSelectors: true },
  checkSourceBySpec: (spec, state, linkerMap) => {
    if (!('kind' in spec) || spec.kind !== 'PColumn') return false;
    const xSource = state.components.x.selectorStates[0]?.selectedSource;
    const ySource = state.components.y.selectorStates[0]?.selectedSource;
    if (!xSource || !ySource) return false;
    const xId = columnOrAxisIdFromString(xSource);
    const yId = columnOrAxisIdFromString(ySource);
    const axisId = isAxisId(xId) ? xId : isAxisId(yId) ? yId : null;   // ← LOAD-BEARING
    if (!axisId) return false;
    const axisSpec = { name: axisId.name, type: axisId.type, domain: axisId.domain } as AxisSpec;
    const normalized = getNormalizedAxesList([axisSpec]);
    const reachable = getNormalizedAxesList(linkerMap.getReachableByLinkersAxesFromAxesNormalized(normalized));
    const allowed = new Set<string>([
      ...normalized.map((a) => AxisId.fromAxisSpec(a).toCanonicalString()),
      ...reachable.map((a) => AxisId.fromAxisSpec(a).toCanonicalString()),
    ]);
    return AxisId.fromAxesSpec(spec.axesSpec).some((id) => allowed.has(id.toCanonicalString()));
  },
}),
```

When neither X nor Y is bound to an axis ID, `axisId` resolves to `null`, the function returns `false` for every candidate, and the dropdown shows nothing. There is no path through this code that examines the axes carried *inside* an X/Y PColumn.

## The semantic fix

What `checkSourceBySpec` should encode:

> *"An additional-curve PColumn must share at least one axis with X's or Y's axis-set, where:*
> - *if X (or Y) is bound to an axis, the axis-set is `{ that axis }`,*
> - *if X (or Y) is bound to a PColumn, the axis-set is the column's `axesSpec`,*
> *plus any axes reachable through the linker map from those starting axes."*

### Recommended implementation

Collect the starting axis specs from *both* X and Y, regardless of whether each is axis-bound or column-bound:

```ts
checkSourceBySpec: (spec, state, linkerMap) => {
  if (!('kind' in spec) || spec.kind !== 'PColumn') return false;
  const xSource = state.components.x.selectorStates[0]?.selectedSource;
  const ySource = state.components.y.selectorStates[0]?.selectedSource;
  if (!xSource || !ySource) return false;

  const startingAxisSpecs: AxisSpec[] = [];
  for (const sourceStr of [xSource, ySource]) {
    const id = columnOrAxisIdFromString(sourceStr);
    if (isAxisId(id)) {
      startingAxisSpecs.push({ name: id.name, type: id.type, domain: id.domain });
    } else if (isColumnId(id)) {
      // NEW path: pull the column's axesSpec via a spec lookup. See
      // "Implementation question" below for how to wire this through.
      const colSpec = linkerMap.getColumnSpecById?.(id);
      if (colSpec?.axesSpec) {
        startingAxisSpecs.push(...colSpec.axesSpec);
      }
    }
  }
  if (startingAxisSpecs.length === 0) return false;

  const normalized = getNormalizedAxesList(startingAxisSpecs);
  const reachable = getNormalizedAxesList(linkerMap.getReachableByLinkersAxesFromAxesNormalized(normalized));
  const allowed = new Set<string>([
    ...normalized.map((a) => AxisId.fromAxisSpec(a).toCanonicalString()),
    ...reachable.map((a) => AxisId.fromAxisSpec(a).toCanonicalString()),
  ]);
  return AxisId.fromAxesSpec(spec.axesSpec).some((id) => allowed.has(id.toCanonicalString()));
},
```

Apply the same change to `scatterplot-umap` if it gains an `additionalCurves` controller too (verify; treat as a parallel fix).

### Implementation question — how does the check fetch a PColumn's `axesSpec`?

The check receives `(spec, state, linkerMap)`. `state` doesn't carry column specs; `linkerMap` already exposes axis-reachability queries, so it's the natural place to add a `getColumnSpecById(id) → PColumnSpec | undefined` accessor.

Investigate the linkerMap factory in `core/visualizations/packages/pf-plots/src/utils/`. `createLinkerMap` (or its modern equivalent in `1.2.0`) almost certainly captures the full set of column specs at construction time — the new accessor would just expose what's already in scope.

If that turns out to require touching too many call sites, the alternative is to thread a separate `specsProvider` through `checkSourceBySpec`'s signature. More invasive — every chart type's `checkSourceBySpec` would need updating. Prefer the linkerMap approach.

## Tests

Add tests in the scatterplot controller's test file (find by grepping `ScatterplotStateController` after the controller exists locally).

| # | X bound to | Y bound to | additionalCurves candidate | Expected | Rationale |
|---|---|---|---|---|---|
| 1 | axis A | column with axes [A, B] | column with axes [A, B] | accept | regression — axis-X case still works |
| 2 | column with axes [A, B] | column with axes [A, B] | column with axes [A, B] | **accept** | new — PColumn-X case the titeseq block needs |
| 3 | column with axes [A, B] | column with axes [A, B] | column with axes [C] | reject | unrelated axis — must still reject |
| 4 | column with axes [A] | column with axes [A, B] | column with axes [A, B] | accept | starting-set is the union of X's and Y's axes |
| 5 | column with axes [A] (where A links to B via linker) | column with axes [A] | column with axes [B] | accept | linker-reachable axis still extends `allowed` |
| 6 | nothing bound | nothing bound | any column | reject | unchanged — early return on missing X/Y stays |

Behavioral, not implementation: assert the return value of `checkSourceBySpec`, not internal helpers.

## Files likely touched (final check against actual layout once `1.2.0` source surfaces)

| File | Change |
|---|---|
| `core/visualizations/packages/pf-plots/src/controllers/controllersByChartType/scatterplot.ts` | Add `additionalCurves` controller (if absent) or fix `checkSourceBySpec` body. |
| `core/visualizations/packages/pf-plots/src/controllers/controllersByChartType/scatterplot-umap.ts` | Same change if the slot exists there too (verify). |
| `core/visualizations/packages/pf-plots/src/utils/linkerMap.ts` (or wherever `createLinkerMap` lives) | Add `getColumnSpecById` accessor to LinkerMap shape. |
| Scatterplot controller test file | Add the six cases above. |
| `core/visualizations/packages/pf-plots/CHANGELOG.md` | Patch-level entry: *"`additionalCurves` source check now accepts PColumn-bound X/Y by examining the column's axesSpec, not just direct axis bindings."* |

## Migration after this PR lands

Once the fix ships in a new `pf-plots` version and titeseq-analysis bumps the catalog:

1. In `worktrees/titeseq-user-feedback-fixes/ui/src/pages/TitrationCurvesPage.vue`, change the X source from the axis spec (the Option A workaround being applied now) back to `concentrationValue.spec`.
2. If the `concentrationValue` PColumn was removed during Option A's cleanup, restore it from git history. The `concentration_value.tsv` Python output, the `concValueColSpec` workflow block, the `--out-concentration-value` CLI arg, and the import wiring all need to come back. Look at the commit immediately preceding the Option A revert as the restore point.
3. Update `docs/investigations/concentration-axis-spec-realignment.md` to add a "Resolved" section noting the SDK fix landed and Option A was reverted.

The audit doc currently flags this fix as "out of scope" for the block's branch. After the fix lands, that section should move to "Resolved by SDK PR #X (pf-plots@<version>)".

## Out of scope

- Other chart types (bubble, dendro, heatmap, etc.) if they have similar PColumn-X gating issues. File separate PRs.
- Adding axis types beyond `Int | Long | String` to the workflow-tengo SDK's axis-type regex — that's a much larger spec-realignment effort, separate from this surgical scatterplot fix.
- Anything in titeseq-analysis itself — handled in the migration step above, *after* this fix lands.

## Risks

- **Linker reachability semantics.** The new code adds *all* of X's and Y's axes to `startingAxisSpecs`, not just one. If `getReachableByLinkersAxesFromAxesNormalized` walks heavily through linkers, the `allowed` set may grow larger than before, accepting some columns the old strict check would have rejected. Probably benign (more permissive, not less) — verify with test case 3.
- **LinkerMap API surface.** Adding `getColumnSpecById` may require updating `createLinkerMap` callers if its return type is narrowly inferred. TypeScript will surface the call sites; update mechanically.
- **Backward compat with snapshots.** No persisted state shape changes — no migration needed for existing project files.
- **The `1.2.0` source mystery.** If the upstream source surfaces during this work, reconcile the two paths: this plan's "feature add against `1.1.66`" approach should stay close to the published `1.2.0` shape so the fix is a straightforward delta, not a rewrite.

## End-to-end validation against titeseq-analysis

After merging the SDK fix and bumping the published version locally:

1. In the titeseq worktree, bump the `@milaboratories/pf-plots` catalog version in `pnpm-workspace.yaml`.
2. `pnpm install && pnpm run build:dev`.
3. Revert Option A in `TitrationCurvesPage.vue` (X source back to `concentrationValue.spec`).
4. `pl mcp update_block` against the TiteSeq desktop project. Run the block.
5. Open Titration Curves. Confirm:
   - X-axis ticks render at user's input concentrations on log scale.
   - **Additional curves** dropdown offers `Fitted mean bin` as a selectable source.
   - Selecting it overlays the Hill sigmoid on each clonotype facet.
   - `validationWarnings: []`, no output errors.

## Quick repro reference for the implementer

To see the bug today (without applying any fix), in the titeseq-analysis branch:

1. `pnpm install && pnpm run build:dev` in `worktrees/titeseq-user-feedback-fixes/`.
2. `pl mcp update_block` against project `NG:0x444c89`, the Tite-Seq Analysis block.
3. Run the block. Open Titration Curves in the desktop app.
4. On the right sidebar, scroll to **Additional curves**. The label is greyed and the dropdown contains no source — even though `fittedMeanBin` exists with the correct shape and is offered in every other slot's picker. **That's the bug.** It's reproducible without any code changes — just a fresh build of the current branch.
