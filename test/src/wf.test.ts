import type {
  BlockData,
  BlockOutputs,
  model,
} from '@platforma-open/platforma-open.titeseq-analysis.model';
import type { InferBlockState, PlRef } from '@platforma-sdk/model';
import { createPlDataTableStateV2, uniquePlId, wrapOutputs } from '@platforma-sdk/model';
import type {
  BlockArgs as SamplesAndDataBlockArgs,
  PlId,
} from '@platforma-open/milaboratories.samples-and-data.model';
import { awaitStableState, blockTest } from '@platforma-sdk/test';
import { blockSpec as myBlockSpec } from 'this-block';
import { blockSpec as samplesAndDataBlockSpec } from '@platforma-open/milaboratories.samples-and-data';
import { blockSpec as importVdjBlockSpec } from '@platforma-open/milaboratories.import-vdj';
import type { BlockArgs as ImportVdjBlockArgs } from '@platforma-open/milaboratories.import-vdj.model';
import { prepareSamplesAndDataArgs } from './fixtures';

blockTest('empty inputs', { timeout: 20000 }, async ({ rawPrj: project, expect }) => {
  const blockId = await project.addBlock('Titeseq Analysis', myBlockSpec);

  const stableState = (await awaitStableState(
    project.getBlockState(blockId),
    15000,
  )) as InferBlockState<typeof model>;

  expect(stableState.outputs).toMatchObject({
    abundanceOptions: { ok: true, value: [] },
    concentrationOptions: { ok: true, value: [] },
    binOptions: { ok: true, value: [] },
    antigenOptions: { ok: true, value: [] },
    validationWarnings: { ok: true, value: [] },
    binMode: { ok: true, value: false },
  });

  const outputs = wrapOutputs<BlockOutputs>(stableState.outputs);
  expect(outputs.abundanceOptions).toEqual([]);
  expect(outputs.concentrationOptions).toEqual([]);
  expect(outputs.binOptions).toEqual([]);
  expect(outputs.antigenOptions).toEqual([]);
  expect(outputs.validationWarnings).toEqual([]);
  expect(outputs.binMode).toEqual(false);

  // Optional workflow outputs appear as slots with value: undefined before
  // the block is runnable (argsValid === false).
  expect(outputs.summaryPfHandle).toBeUndefined();
  expect(outputs.titrationCurvesPf).toBeUndefined();
  expect(outputs.logHandle).toBeUndefined();

  // isRunning should be falsy (undefined or false) while the block is idle.
  expect(outputs.isRunning).toBeFalsy();
});

blockTest(
  'bin_mode fixture: option lists + end-to-end summary contract',
  { timeout: 900000 },
  async ({ rawPrj: project, ml, helpers, expect }) => {
    // One S&D ingest powers both the option-list checks and the end-to-end run
    // below — each blockTest spins up a fresh platforma container, so merging
    // these tests cuts fixture uploads roughly in half (CI runner disk is the
    // bottleneck; see samples-and-data/test for the sample-count pattern).
    const sndBlockId = await project.addBlock('Samples & Data', samplesAndDataBlockSpec);
    const importBlockId = await project.addBlock('Import V(D)J Data', importVdjBlockSpec);
    const titeseqBlockId = await project.addBlock('Titeseq Analysis', myBlockSpec);

    const { args: sndArgs, loaded } = await prepareSamplesAndDataArgs('bin_mode', helpers);
    await project.setBlockArgs(sndBlockId, sndArgs);
    await project.runBlock(sndBlockId);
    await helpers.awaitBlockDone(sndBlockId, 60000);

    // import-vdj-data exposes the S&D dataset via its `datasetOptions` retentive
    // output once S&D is done. Pick the first (only) option for this fixture.
    const importState1 = await awaitStableState(
      project.getBlockState(importBlockId),
      30000,
    );
    expect(importState1.outputs).toBeDefined();
    const importOutputs1 = wrapOutputs(importState1.outputs!) as {
      datasetOptions?: { ref: PlRef; label: string }[];
    };
    expect(importOutputs1.datasetOptions).toBeDefined();
    expect(importOutputs1.datasetOptions!.length).toBeGreaterThan(0);
    const datasetRef = importOutputs1.datasetOptions![0].ref;

    await project.setBlockArgs(importBlockId, {
      defaultBlockLabel: '',
      customBlockLabel: '',
      datasetRef,
      format: 'mixcr',
      chains: ['IGHeavy'],
    } satisfies ImportVdjBlockArgs);
    await project.runBlock(importBlockId);
    await helpers.awaitBlockDone(importBlockId, 300000);

    // Tier-1: option lists populate correctly once import-vdj-data emits the
    // isAnchor abundance column and the per-sample metadata columns propagate.
    const idleState = (await awaitStableState(
      project.getBlockState(titeseqBlockId),
      30000,
    )) as InferBlockState<typeof model>;
    const idleOutputs = wrapOutputs<BlockOutputs>(idleState.outputs);

    // Abundance: exactly one isAnchor [sampleId, clonotypeKey] PColumn from
    // import-vdj-data. A second entry here would mean the filter is matching
    // an unintended column.
    expect(idleOutputs.abundanceOptions.length).toBe(1);

    // Concentration/bin options also include per-sample stats PColumns emitted
    // by import-vdj-data (clones count, total reads, etc.) — assert the
    // S&D-sourced metadata column is present rather than pinning the count.
    const concLabels = idleOutputs.concentrationOptions.map((o) => o.label);
    const binLabels = idleOutputs.binOptions.map((o) => o.label);
    expect(concLabels).toEqual(expect.arrayContaining(['antigen_conc_M']));
    expect(binLabels).toEqual(expect.arrayContaining(['bin_number']));
    // No antigen metadata in the bin_mode fixture — the fixture's antigen
    // label must not appear. (S&D emits a `Sample Name` String column that
    // matches the antigen filter shape, so we can't assert the list is empty.)
    const antigenLabels = idleOutputs.antigenOptions.map((o) => o.label);
    expect(antigenLabels).not.toContain('antigen');

    // Filter discrimination: the integer bin column must not leak into the
    // concentration (Float/Double) dropdown, and vice versa.
    const binRefs = new Set(idleOutputs.binOptions.map((o) => JSON.stringify(o.ref)));
    for (const concOption of idleOutputs.concentrationOptions) {
      expect(binRefs.has(JSON.stringify(concOption.ref))).toBe(false);
    }

    // Titeseq block is idle before we configure it; workflow outputs must not
    // appear yet. Cross-reference the metadata ids against the fixture helper
    // so a plumbing mismatch surfaces before we spend time on the run.
    expect(idleOutputs.validationWarnings).toEqual([]);
    expect(idleOutputs.isRunning).toBeFalsy();
    expect(idleOutputs.summaryPfHandle).toBeUndefined();
    expect(idleOutputs.titrationCurvesPf).toBeUndefined();
    expect(Object.keys(loaded.sampleIdByName).length).toBeGreaterThan(0);
    expect(loaded.metadataIds.concentration).toBeDefined();
    expect(loaded.metadataIds.bin).toBeDefined();

    const abundanceRef = idleOutputs.abundanceOptions[0].ref;
    const concentrationRef = idleOutputs.concentrationOptions.find(
      (o) => o.label === 'antigen_conc_M',
    )!.ref;
    const binRef = idleOutputs.binOptions.find((o) => o.label === 'bin_number')!.ref;

    // Configure the titeseq block with bin-mode args. DataModelBuilder merges
    // args + UI state into a single BlockData record (see model/src/index.ts);
    // mutateBlockStorage with `update-block-data` is the V3 equivalent of
    // setBlockArgs + setUiState.
    await project.mutateBlockStorage(titeseqBlockId, {
      operation: 'update-block-data',
      value: {
        abundanceRef,
        concentrationColumnRef: concentrationRef,
        binColumnRef: binRef,
        antigenColumnRef: undefined,
        targetAntigen: undefined,
        minReadsPerConcentration: 3,
        minConcentrationPoints: 5,
        r2ThresholdGood: 0.8,
        r2ThresholdFailed: 0.5,
        nMin: 0.5,
        nMax: 2.0,
        hookEffectThresholdBin: 0.2,
        hookEffectThresholdNoBin: 0.02,
        hookEffectMinReads: 20,
        defaultBlockLabel: 'Tite-Seq Analysis',
        customBlockLabel: '',
        tableState: createPlDataTableStateV2(),
        graphStateTitrationCurves: {
          title: 'Titration Curves',
          template: 'dots',
          currentTab: null,
          layersSettings: { dots: {} },
          axesSettings: { axisX: { scale: 'log' } },
        },
        graphStateKDHistogram: {
          title: 'Kd,app Distribution',
          template: 'bins',
          currentTab: null,
          layersSettings: { bins: {} },
          axesSettings: { axisX: { scale: 'log' }, other: { binsCount: 30 } },
        },
        graphStateAffinityVsFit: {
          title: 'Affinity vs Fit Quality',
          template: 'dots',
          currentTab: null,
          layersSettings: { dots: {} },
          axesSettings: { axisX: { scale: 'log' } },
        },
        settingsOpen: false,
      } satisfies BlockData,
    });

    await project.runBlock(titeseqBlockId);
    const doneState = await helpers.awaitBlockDoneAndGetStableBlockState(
      titeseqBlockId,
      360000,
    );
    const outputs = wrapOutputs<BlockOutputs>(
      doneState.outputs as unknown as BlockOutputs,
    );

    // Tier-3 contract: summary PColumn specs must all be present and match the
    // expected names/valueTypes. Downstream consumers (Lead Selection, UI) key
    // off these exact names, so any rename must be deliberate.
    expect(outputs.summaryPfCols).toBeDefined();
    const specByName = new Map(outputs.summaryPfCols!.map((c) => [c.spec.name, c.spec]));
    const expectedCols: Array<[string, string]> = [
      ['pl7.app/vdj/kd', 'Double'],
      ['pl7.app/vdj/hillCoefficient', 'Double'],
      ['pl7.app/vdj/curveFitR2', 'Double'],
      ['pl7.app/vdj/affinityClass', 'String'],
      ['pl7.app/vdj/fitFailureReason', 'String'],
      ['pl7.app/vdj/kdOutOfRange', 'String'],
    ];
    for (const [name, valueType] of expectedCols) {
      const spec = specByName.get(name);
      expect(spec, `summary column ${name}`).toBeDefined();
      expect(spec!.valueType, `${name} valueType`).toBe(valueType);
    }

    // Every summary column must be axed by a single clonotypeKey axis.
    for (const [name, spec] of specByName) {
      expect(spec.axesSpec.length, `${name} axes count`).toBe(1);
      expect(spec.axesSpec[0].name, `${name} axis name`).toBe('pl7.app/vdj/clonotypeKey');
    }

    // Tier-2 happy path: the summary PFrame must contain data, and the
    // affinityClass distribution must roughly track the manifest (bin_mode has
    // 5 clonotypes → 2 Good, 1 Partial, 2 Failed). Lenient comparison — the
    // Hill fitter has some tolerance baked in, so we assert "at least one Good"
    // rather than pinning exact counts. Data extraction goes through the
    // `summaryTable` PTable handle; PFrame handles don't support getSpec/getData.
    expect(outputs.summaryPfHandle).toBeDefined();
    const summaryTable = outputs.summaryTable;
    if (!summaryTable) throw new Error('summaryTable unexpectedly undefined after block ran');
    const tableHandle = summaryTable.fullTableHandle;
    if (!tableHandle) throw new Error('summaryTable.fullTableHandle unexpectedly undefined after block ran');
    const tableSpec = await ml.driverKit.pFrameDriver.getSpec(tableHandle);
    const affinityIdx = tableSpec.findIndex(
      (c) => c.type === 'column' && c.spec.name === 'pl7.app/vdj/affinityClass',
    );
    const kdIdx = tableSpec.findIndex(
      (c) => c.type === 'column' && c.spec.name === 'pl7.app/vdj/kd',
    );
    const r2Idx = tableSpec.findIndex(
      (c) => c.type === 'column' && c.spec.name === 'pl7.app/vdj/curveFitR2',
    );
    expect(affinityIdx, 'affinityClass column').toBeGreaterThanOrEqual(0);
    expect(kdIdx, 'kd column').toBeGreaterThanOrEqual(0);
    expect(r2Idx, 'curveFitR2 column').toBeGreaterThanOrEqual(0);

    const tableData = await ml.driverKit.pFrameDriver.getData(
      tableHandle,
      [affinityIdx, kdIdx, r2Idx],
    );
    const affinityValues = tableData[0].data as unknown as (string | null | undefined)[];
    const kdValues = tableData[1].data as unknown as (number | null | undefined)[];

    expect(affinityValues.length, 'clonotype row count').toBeGreaterThanOrEqual(5);

    const affinityCounts = { Good: 0, Partial: 0, Failed: 0, other: 0 };
    for (const v of affinityValues) {
      if (v === 'Good') affinityCounts.Good += 1;
      else if (v === 'Partial') affinityCounts.Partial += 1;
      else if (v === 'Failed') affinityCounts.Failed += 1;
      else affinityCounts.other += 1;
    }
    expect(affinityCounts.Good, `Good count (distribution: ${JSON.stringify(affinityCounts)})`).toBeGreaterThanOrEqual(1);
    expect(affinityCounts.Failed, `Failed count (distribution: ${JSON.stringify(affinityCounts)})`).toBeGreaterThanOrEqual(1);

    // Tier-2 numeric check: at least one Kd in the expected range for Good
    // clonotypes. G_LOW targets [3e-10, 3e-09], G_MID targets [3e-09, 3e-08],
    // so any finite Kd in [3e-10, 3e-08] is consistent with a Good fit.
    const goodKds: number[] = [];
    for (let i = 0; i < affinityValues.length; i += 1) {
      if (affinityValues[i] === 'Good') {
        const kd = kdValues[i];
        if (typeof kd === 'number' && Number.isFinite(kd)) goodKds.push(kd);
      }
    }
    expect(goodKds.length, 'Good clonotypes with finite Kd').toBeGreaterThanOrEqual(1);
    const inRange = goodKds.filter((kd) => kd >= 3e-10 && kd <= 3e-08);
    expect(inRange.length, `Kd in expected range (Good Kds: ${goodKds.join(', ')})`).toBeGreaterThanOrEqual(1);

    // titrationCurvesPf combines summary + signal columns, so it must carry
    // all 6 summary names plus meanBin + fittedMeanBin. Use listColumns here
    // because titrationCurvesPf is a PFrame (no getSpec/getData).
    const titrationCurvesPf = outputs.titrationCurvesPf;
    if (!titrationCurvesPf) throw new Error('titrationCurvesPf unexpectedly undefined after block ran');
    const curvesCols = await ml.driverKit.pFrameDriver.listColumns(titrationCurvesPf);
    const curveNames = new Set(curvesCols.map((c) => c.spec.name));
    for (const [name] of expectedCols) {
      expect(curveNames.has(name), `titrationCurvesPf missing ${name}`).toBe(true);
    }
    expect(curveNames.has('pl7.app/vdj/meanBin'), 'meanBin').toBe(true);
    expect(curveNames.has('pl7.app/vdj/fittedMeanBin'), 'fittedMeanBin').toBe(true);

    // isEmpty should be false — we produced at least 5 rows.
    expect(outputs.isEmpty).toBe(false);
    // After the run settles, the block is idle again.
    expect(outputs.isRunning).toBeFalsy();
  },
);

blockTest(
  'bin_mode + sort_fraction: FACS correction plumbs through and toggles annotation',
  { timeout: 900000 },
  async ({ rawPrj: project, helpers, expect }) => {
    // Goal: exercise the full sort_fraction wiring (model → workflow → Python →
    // signal annotation → model output → UI badge). Correctness is covered by
    // Python unit/integration/regression tests; this test is a plumbing guard.
    // Uniform 0.25 per bin is a mathematical no-op but still verifies the
    // correction code path runs end-to-end and the annotation round-trips.
    const sndBlockId = await project.addBlock('Samples & Data', samplesAndDataBlockSpec);
    const importBlockId = await project.addBlock('Import V(D)J Data', importVdjBlockSpec);
    const titeseqBlockId = await project.addBlock('Titeseq Analysis', myBlockSpec);

    const { args: baseSndArgs, loaded } = await prepareSamplesAndDataArgs('bin_mode', helpers);

    // Inject a sort_fraction metadata column: 4 bins per concentration, uniform
    // 0.25 — validator accepts (sums to 1.0 per concentration, in [0,1]).
    const sortFractionId = uniquePlId() as unknown as PlId;
    const sortFractionData: Record<PlId, number> = {};
    for (const name of loaded.raw.sample_ids) {
      sortFractionData[loaded.sampleIdByName[name]] = 0.25;
    }
    const sndArgs: SamplesAndDataBlockArgs = {
      ...baseSndArgs,
      metadata: [
        ...baseSndArgs.metadata,
        {
          id: sortFractionId,
          label: 'sort_fraction',
          global: true,
          valueType: 'Double',
          data: sortFractionData,
        },
      ],
    };

    await project.setBlockArgs(sndBlockId, sndArgs);
    await project.runBlock(sndBlockId);
    await helpers.awaitBlockDone(sndBlockId, 60000);

    const importState = await awaitStableState(project.getBlockState(importBlockId), 30000);
    expect(importState.outputs).toBeDefined();
    const importOutputs = wrapOutputs(importState.outputs!) as {
      datasetOptions?: { ref: PlRef; label: string }[];
    };
    const datasetRef = importOutputs.datasetOptions![0].ref;
    await project.setBlockArgs(importBlockId, {
      defaultBlockLabel: '',
      customBlockLabel: '',
      datasetRef,
      format: 'mixcr',
      chains: ['IGHeavy'],
    } satisfies ImportVdjBlockArgs);
    await project.runBlock(importBlockId);
    await helpers.awaitBlockDone(importBlockId, 300000);

    const idleState = (await awaitStableState(
      project.getBlockState(titeseqBlockId),
      30000,
    )) as InferBlockState<typeof model>;
    const idleOutputs = wrapOutputs<BlockOutputs>(idleState.outputs);

    // The sort_fraction column must surface in the dedicated dropdown (Float
    // [sampleId] filter) — not in the bin dropdown (Long).
    const sortFractionLabels = idleOutputs.sortFractionOptions?.map((o) => o.label) ?? [];
    expect(sortFractionLabels).toEqual(expect.arrayContaining(['sort_fraction']));

    const abundanceRef = idleOutputs.abundanceOptions[0].ref;
    const concentrationRef = idleOutputs.concentrationOptions.find(
      (o) => o.label === 'antigen_conc_M',
    )!.ref;
    const binRef = idleOutputs.binOptions.find((o) => o.label === 'bin_number')!.ref;
    const sortFractionRef = idleOutputs.sortFractionOptions!.find(
      (o) => o.label === 'sort_fraction',
    )!.ref;

    await project.mutateBlockStorage(titeseqBlockId, {
      operation: 'update-block-data',
      value: {
        abundanceRef,
        concentrationColumnRef: concentrationRef,
        binColumnRef: binRef,
        sortFractionColumnRef: sortFractionRef,
        antigenColumnRef: undefined,
        targetAntigen: undefined,
        minReadsPerConcentration: 3,
        minConcentrationPoints: 5,
        r2ThresholdGood: 0.8,
        r2ThresholdFailed: 0.5,
        nMin: 0.5,
        nMax: 2.0,
        hookEffectThresholdBin: 0.2,
        hookEffectThresholdNoBin: 0.02,
        hookEffectMinReads: 20,
        defaultBlockLabel: 'Tite-Seq Analysis',
        customBlockLabel: '',
        tableState: createPlDataTableStateV2(),
        graphStateTitrationCurves: {
          title: 'Titration Curves',
          template: 'dots',
          currentTab: null,
          layersSettings: { dots: {} },
          axesSettings: { axisX: { scale: 'log' } },
        },
        graphStateKDHistogram: {
          title: 'Kd,app Distribution',
          template: 'bins',
          currentTab: null,
          layersSettings: { bins: {} },
          axesSettings: { axisX: { scale: 'log' }, other: { binsCount: 30 } },
        },
        graphStateAffinityVsFit: {
          title: 'Affinity vs Fit Quality',
          template: 'dots',
          currentTab: null,
          layersSettings: { dots: {} },
          axesSettings: { axisX: { scale: 'log' } },
        },
        settingsOpen: false,
      } satisfies BlockData,
    });

    await project.runBlock(titeseqBlockId);
    const doneState = await helpers.awaitBlockDoneAndGetStableBlockState(
      titeseqBlockId,
      360000,
    );
    const outputs = wrapOutputs<BlockOutputs>(doneState.outputs as unknown as BlockOutputs);

    // Plumbing assertions per plan §12e: the annotation round-trips from
    // workflow/src/main.tpl.tengo through the signalPf → facsCorrectionActive
    // output → UI badge.
    expect(outputs.facsCorrectionActive).toBe(true);

    // meanBin PColumn must be present in titrationCurvesPf (already asserted in
    // the plain bin_mode test) AND carry the facsCorrected annotation on the
    // signalPf side that `facsCorrectionActive` reads from.
    expect(outputs.titrationCurvesPf).toBeDefined();

    // Block completed cleanly — no validation warnings from the sort_fraction
    // validator. (The validator runs server-side in Python; any malformed
    // sort-fraction setup would surface as a block error, not a warning here.)
    expect(outputs.validationWarnings).toEqual([]);
    expect(outputs.isEmpty).toBe(false);
    expect(outputs.isRunning).toBeFalsy();
  },
);
