// Reproductions from the 2026-10-20 bug hunt. Each test asserts the correct
// behaviour; a test marked `.fails` pins a defect that is still open.
import type {
  BlockData,
  BlockOutputs,
  platforma,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import type { InferBlockState, PlRef } from "@platforma-sdk/model";
import { createPlDataTableStateV2, uniquePlId, wrapOutputs } from "@platforma-sdk/model";
import type {
  BlockArgs as SamplesAndDataBlockArgs,
  PlId,
} from "@platforma-open/milaboratories.samples-and-data.model";
import { awaitStableState, blockTest } from "@platforma-sdk/test";
import { TiteseqAnalysisBlockPointer as myBlockSpec } from "this-block";
import { blockSpec as samplesAndDataBlockSpec } from "@platforma-open/milaboratories.samples-and-data";
import { blockSpec as importVdjBlockSpec } from "@platforma-open/milaboratories.import-vdj";
import type { BlockArgs as ImportVdjBlockArgs } from "@platforma-open/milaboratories.import-vdj.model";
import { prepareSamplesAndDataArgs } from "./fixtures";

type TestContext = Parameters<NonNullable<Parameters<typeof blockTest>[2]>>[0];
type Project = TestContext["rawPrj"];
type Helpers = TestContext["helpers"];

function blockData(over: Partial<BlockData>): BlockData {
  return {
    abundanceRef: undefined,
    concentrationColumnRef: undefined,
    binColumnRef: undefined,
    antigenColumnRef: undefined,
    sortFractionColumnRef: undefined,
    cellCountColumnRef: undefined,
    targetAntigen: undefined,
    // The fixtures were built on the linear Hill model.
    fitModelOverride: "hill_free_n",
    binFluorescenceColumnRef: undefined,
    gateWidthColumnRef: undefined,
    amplitudePriorQ: undefined,
    // The fixtures' expected classes come from a free-n Hill fit with the hook detector on,
    // which is no longer the default.
    hookDetectorEnabled: true,
    pinHillExponent: false,
    binWidthDecades: 0.65,
    minPlateauFraction: 0.65,
    baselineFluorescence: undefined,
    kdSdGood: 0.5,
    minReadsPerConcentration: 3,
    minConcentrationPoints: 5,
    r2ThresholdGood: 0.8,
    r2ThresholdFailed: 0.5,
    nMin: 0.5,
    nMax: 2.0,
    hookEffectThresholdBin: 0.2,
    hookEffectThresholdNoBin: 0.02,
    hookEffectMinReads: 20,
    defaultBlockLabel: "Tite-Seq Analysis",
    customBlockLabel: "",
    tableState: createPlDataTableStateV2(),
    graphStateTitrationCurves: {
      title: "",
      template: "dots",
      currentTab: null,
      layersSettings: { dots: {} },
      axesSettings: { axisX: { scale: "log" } },
    },
    graphStateKDHistogram: {
      title: "",
      template: "bins",
      currentTab: null,
      layersSettings: { bins: {} },
      axesSettings: { axisX: { scale: "log" }, other: { binsCount: 30 } },
    },
    graphStateAffinityVsFit: {
      title: "",
      template: "dots",
      currentTab: null,
      layersSettings: { dots: {} },
      axesSettings: { axisX: { scale: "log" } },
    },
    settingsOpen: false,
    ...over,
  };
}

/**
 * S&D (antigen fixture, plus a 0.25 sort fraction and an integer pM
 * concentration column) -> Import V(D)J -> Titeseq. Returns the refs the
 * Titeseq block offers once the upstream blocks are done.
 */
async function setupAntigenChain(project: Project, helpers: Helpers) {
  const sndBlockId = await project.addBlock("Samples & Data", samplesAndDataBlockSpec);
  const importBlockId = await project.addBlock("Import V(D)J Data", importVdjBlockSpec);
  const titeseqBlockId = await project.addBlock("Titeseq Analysis", myBlockSpec);

  const { args: baseSndArgs, loaded } = await prepareSamplesAndDataArgs("antigen", helpers);
  const sortFraction: Record<PlId, number> = {};
  const concPm: Record<PlId, number> = {};
  for (const name of loaded.raw.sample_ids) {
    const id = loaded.sampleIdByName[name];
    sortFraction[id] = 0.25;
    concPm[id] = Math.round(loaded.raw.concentration[name] * 1e12);
  }
  const sndArgs: SamplesAndDataBlockArgs = {
    ...baseSndArgs,
    metadata: [
      ...baseSndArgs.metadata,
      {
        id: uniquePlId() as unknown as PlId,
        label: "sort_fraction",
        global: true,
        valueType: "Double",
        data: sortFraction,
      },
      {
        id: uniquePlId() as unknown as PlId,
        label: "conc_pM",
        global: true,
        valueType: "Long",
        data: concPm,
      },
    ],
  };
  await project.setBlockArgs(sndBlockId, sndArgs);
  await project.runBlock(sndBlockId);
  await helpers.awaitBlockDone(sndBlockId, 120000);

  const importState = await awaitStableState(project.getBlockState(importBlockId), 30000);
  const importOutputs = wrapOutputs(importState.outputs!) as {
    datasetOptions?: { ref: PlRef; label: string }[];
  };
  await project.setBlockArgs(importBlockId, {
    defaultBlockLabel: "",
    customBlockLabel: "",
    datasetRef: importOutputs.datasetOptions![0].ref,
    format: "mixcr",
    chains: ["IGHeavy"],
  } satisfies ImportVdjBlockArgs);
  await project.runBlock(importBlockId);
  await helpers.awaitBlockDone(importBlockId, 300000);

  const idle = wrapOutputs<BlockOutputs>(
    (
      (await awaitStableState(project.getBlockState(titeseqBlockId), 30000)) as InferBlockState<
        typeof platforma
      >
    ).outputs,
  );
  const byLabel = (opts: { ref: PlRef; label: string }[], label: string) => {
    const hit = opts.find((o) => o.label === label);
    if (!hit) throw new Error(`setup: option "${label}" not offered`);
    return hit.ref;
  };
  return {
    titeseqBlockId,
    idle,
    refs: {
      abundance: idle.abundanceOptions[0].ref,
      concentration: byLabel(idle.concentrationOptions, "antigen_conc_M"),
      bin: byLabel(idle.binOptions, "bin_number"),
      sortFraction: byLabel(idle.sortFractionOptions ?? [], "sort_fraction"),
      antigen: byLabel(idle.antigenOptions, "antigen"),
    },
  };
}

blockTest(
  "bug-hunt H4: sort fraction plus antigen column fails the run (sums checked across antigens)",
  { timeout: 900000 },
  async ({ rawPrj: project, helpers, expect }) => {
    const { titeseqBlockId, refs } = await setupAntigenChain(project, helpers);
    await project.mutateBlockStorage(titeseqBlockId, {
      operation: "update-block-data",
      value: blockData({
        abundanceRef: refs.abundance,
        concentrationColumnRef: refs.concentration,
        binColumnRef: refs.bin,
        antigenColumnRef: refs.antigen,
        targetAntigen: "target",
        sortFractionColumnRef: refs.sortFraction,
      }),
    });
    await project.runBlock(titeseqBlockId);
    // Every antigen's bins sum to 1.0 per concentration; the block must accept it.
    await expect(helpers.awaitBlockDone(titeseqBlockId, 360000)).resolves.toBeUndefined();
  },
);

blockTest(
  "bug-hunt M5: validationWarnings itself errors when the run fails, so the in-block alert never shows",
  { timeout: 900000 },
  async ({ rawPrj: project, helpers, expect }) => {
    const { titeseqBlockId, refs } = await setupAntigenChain(project, helpers);
    await project.mutateBlockStorage(titeseqBlockId, {
      operation: "update-block-data",
      value: blockData({
        abundanceRef: refs.abundance,
        concentrationColumnRef: refs.concentration,
        binColumnRef: refs.bin,
        antigenColumnRef: refs.antigen,
        targetAntigen: "target",
        sortFractionColumnRef: refs.sortFraction,
        // A 0 M baseline above the plateau prior: the Adams fit rejects the run in Python.
        fitModelOverride: "poisson_log",
        baselineFluorescence: 1000,
      }),
    });
    await project.runBlock(titeseqBlockId);
    await helpers.awaitBlockDone(titeseqBlockId, { timeout: 360000, ignoreBlockError: true });
    const state = (await awaitStableState(
      project.getBlockState(titeseqBlockId),
      30000,
    )) as InferBlockState<typeof platforma>;
    const warnings = state.outputs?.validationWarnings;
    expect(warnings?.ok).toBe(true);
    const messages = warnings?.ok ? warnings.value.map((w) => w.message) : [];
    expect(messages.some((m) => m.includes("plateau"))).toBe(true);
  },
);

blockTest(
  "bug-hunt M6: isEmpty stays false when every clonotype is Failed, so the 'All clonotypes failed' alert never shows",
  { timeout: 900000 },
  async ({ rawPrj: project, helpers, expect }) => {
    const { titeseqBlockId, refs } = await setupAntigenChain(project, helpers);
    await project.mutateBlockStorage(titeseqBlockId, {
      operation: "update-block-data",
      value: blockData({
        abundanceRef: refs.abundance,
        concentrationColumnRef: refs.concentration,
        binColumnRef: refs.bin,
        antigenColumnRef: refs.antigen,
        targetAntigen: "target",
        // Both thresholds at 1.0: no fit reaches R² = 1, so every clonotype is Failed.
        r2ThresholdGood: 1,
        r2ThresholdFailed: 1,
      }),
    });
    await project.runBlock(titeseqBlockId);
    const done = await helpers.awaitBlockDoneAndGetStableBlockState(titeseqBlockId, 360000);
    const outputs = wrapOutputs<BlockOutputs>(done.outputs as unknown as BlockOutputs);
    expect(outputs.isEmpty).toBe(true);
  },
);

blockTest(
  "bug-hunt M7: binMode follows unsaved settings, not the results on screen",
  { timeout: 900000 },
  async ({ rawPrj: project, helpers, expect }) => {
    const { titeseqBlockId, refs } = await setupAntigenChain(project, helpers);
    const ran = {
      abundanceRef: refs.abundance,
      concentrationColumnRef: refs.concentration,
      binColumnRef: refs.bin,
      antigenColumnRef: refs.antigen,
      targetAntigen: "target",
    };
    await project.mutateBlockStorage(titeseqBlockId, {
      operation: "update-block-data",
      value: blockData(ran),
    });
    await project.runBlock(titeseqBlockId);
    await helpers.awaitBlockDone(titeseqBlockId, 360000);

    // Clear the bin column without re-running: the bin-mode results stay on screen.
    await project.mutateBlockStorage(titeseqBlockId, {
      operation: "update-block-data",
      value: blockData({ ...ran, binColumnRef: undefined }),
    });
    const state = wrapOutputs<BlockOutputs>(
      (
        (await awaitStableState(project.getBlockState(titeseqBlockId), 30000)) as InferBlockState<
          typeof platforma
        >
      ).outputs,
    );
    expect(state.summaryPfCols).toBeDefined();
    expect(state.binMode).toBe(true);
  },
);

// Open: concentration options take Float/Double only, which is also what keeps the integer FACS
// bin column out of that list (see the guard in wf.test.ts). Accepting Int/Long fixes this but
// lists bin and read-count columns as concentrations; needs a narrower rule.
blockTest.fails(
  "bug-hunt M4: an integer (Long) concentration column is not offered as a concentration",
  { timeout: 900000 },
  async ({ rawPrj: project, helpers, expect }) => {
    const { idle } = await setupAntigenChain(project, helpers);
    expect(idle.concentrationOptions.map((o) => o.label)).toContain("conc_pM");
  },
);
