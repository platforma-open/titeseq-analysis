import type { GraphMakerState } from "@milaboratories/graph-maker";
import type {
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlDataTableStateV2,
  PlRef,
} from "@platforma-sdk/model";
import {
  BlockModelV3,
  DataModelBuilder,
  createPFrameForGraphs,
  createPlDataTableStateV2,
  createPlDataTableV2,
  isPColumnSpec,
} from "@platforma-sdk/model";

export type ValidationIssue = {
  severity: "warning" | "error";
  message: string;
};

export type BlockArgs = {
  // Inputs (R1–R4)
  abundanceRef?: PlRef;
  concentrationColumnRef?: PlRef;
  binColumnRef?: PlRef;
  antigenColumnRef?: PlRef;
  targetAntigen?: string;

  // Fitting params
  minReadsPerConcentration: number;
  minConcentrationPoints: number;
  r2ThresholdGood: number;
  r2ThresholdFailed: number;
  nMin: number;
  nMax: number;
  hookEffectThresholdBin: number;
  hookEffectThresholdNoBin: number;
  hookEffectMinReads: number;

  // Label
  defaultBlockLabel: string;
  customBlockLabel: string;
};

export type BlockData = BlockArgs & {
  tableState: PlDataTableStateV2;
  graphStateTitrationCurves: GraphMakerState;
  graphStateKDHistogram: GraphMakerState;
  graphStateAffinityVsFit: GraphMakerState;
  settingsOpen: boolean;
};

type LegacyUiState = {
  tableState: PlDataTableStateV2;
  graphStateTitrationCurves: GraphMakerState;
  graphStateKDHistogram: GraphMakerState;
  graphStateAffinityVsFit: GraphMakerState;
  settingsOpen: boolean;
};

const dataModel = new DataModelBuilder()
  .from<BlockData>("v1")
  .upgradeLegacy<BlockArgs, LegacyUiState>(({ args, uiState }) => ({
    ...args,
    ...uiState,
  }))
  .init(() => ({
    abundanceRef: undefined,
    concentrationColumnRef: undefined,
    binColumnRef: undefined,
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
    defaultBlockLabel: "Titeseq Analysis",
    customBlockLabel: "",
    tableState: createPlDataTableStateV2(),
    graphStateTitrationCurves: {
      title: "Titration Curves",
      template: "dots",
      currentTab: null,
      layersSettings: {
        dots: {},
      },
      axesSettings: {
        axisX: {
          scale: "log",
        },
      },
    },
    graphStateKDHistogram: {
      title: "K_D,app Distribution",
      template: "bins",
      currentTab: null,
      layersSettings: {
        bins: {},
      },
      axesSettings: {
        axisX: {
          scale: "log",
        },
        other: { binsCount: 30 },
      },
    },
    graphStateAffinityVsFit: {
      title: "Affinity vs Fit Quality",
      template: "dots",
      currentTab: null,
      layersSettings: {
        dots: {},
      },
      axesSettings: {
        axisX: {
          scale: "log",
        },
      },
    },
    settingsOpen: false,
  }));

function isIntegerValueType(vt: string | undefined): boolean {
  return vt === "Int" || vt === "Long";
}

function isFloatValueType(vt: string | undefined): boolean {
  return vt === "Float" || vt === "Double";
}

export const model = BlockModelV3.create(dataModel)

  .args<BlockArgs>((data) => {
    if (data.abundanceRef === undefined) throw new Error("Abundance column is required");
    if (data.concentrationColumnRef === undefined)
      throw new Error("Concentration column is required");
    if (data.antigenColumnRef !== undefined && !data.targetAntigen)
      throw new Error("Target antigen is required when an antigen column is selected");
    if (data.r2ThresholdFailed > data.r2ThresholdGood)
      throw new Error("Failed R² threshold must be ≤ Good R² threshold");
    if (data.nMin >= data.nMax) throw new Error("Hill coefficient nMin must be < nMax");
    if (data.r2ThresholdGood < 0 || data.r2ThresholdGood > 1)
      throw new Error("R² threshold (Good) out of range");
    if (data.r2ThresholdFailed < 0 || data.r2ThresholdFailed > 1)
      throw new Error("R² threshold (Failed) out of range");
    if (data.minReadsPerConcentration < 1)
      throw new Error("Min reads per concentration must be ≥ 1");
    if (data.minConcentrationPoints < 3) throw new Error("Min concentration points must be ≥ 3");
    return {
      abundanceRef: data.abundanceRef,
      concentrationColumnRef: data.concentrationColumnRef,
      binColumnRef: data.binColumnRef,
      antigenColumnRef: data.antigenColumnRef,
      targetAntigen: data.targetAntigen,
      minReadsPerConcentration: data.minReadsPerConcentration,
      minConcentrationPoints: data.minConcentrationPoints,
      r2ThresholdGood: data.r2ThresholdGood,
      r2ThresholdFailed: data.r2ThresholdFailed,
      nMin: data.nMin,
      nMax: data.nMax,
      hookEffectThresholdBin: data.hookEffectThresholdBin,
      hookEffectThresholdNoBin: data.hookEffectThresholdNoBin,
      hookEffectMinReads: data.hookEffectMinReads,
      defaultBlockLabel: data.defaultBlockLabel,
      customBlockLabel: data.customBlockLabel,
    };
  })

  .output("abundanceOptions", (ctx) =>
    ctx.resultPool.getOptions(
      [
        {
          axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/vdj/clonotypeKey" }],
          annotations: { "pl7.app/isAnchor": "true" },
        },
        {
          axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/vdj/scClonotypeKey" }],
          annotations: { "pl7.app/isAnchor": "true" },
        },
      ],
      {
        label: { includeNativeLabel: false },
      },
    ),
  )

  .output("concentrationOptions", (ctx) => {
    const candidates = ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        isFloatValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    );
    // Hide columns whose values are entirely null — selecting one throws an
    // "empty column" error inside the Python pipeline, so the option should
    // not appear in the dropdown in the first place. Data-not-ready is shown
    // conservatively to avoid options flickering in/out while upstream loads.
    return candidates.filter((opt) => {
      const data = ctx.resultPool.getDataByRef(opt.ref)?.data;
      if (!data) return true;
      const values = data.getDataAsJson<Record<string, number | null>>()?.["data"];
      if (!values) return true;
      return Object.values(values).some((v) => v !== null && v !== undefined);
    });
  })

  .output("binOptions", (ctx) => {
    const candidates = ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        isIntegerValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    );
    // Same reasoning as concentrationOptions: an all-null integer column
    // arrives at the Python side as a String column of empty strings and
    // fails validation, so hide it here.
    return candidates.filter((opt) => {
      const data = ctx.resultPool.getDataByRef(opt.ref)?.data;
      if (!data) return true;
      const values = data.getDataAsJson<Record<string, number | null>>()?.["data"];
      if (!values) return true;
      return Object.values(values).some((v) => v !== null && v !== undefined);
    });
  })

  .output("antigenOptions", (ctx) =>
    ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        spec.valueType === "String" &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    ),
  )

  .output("targetAntigenValues", (ctx): string[] | undefined => {
    if (!ctx.data.antigenColumnRef) return undefined;
    const data = ctx.resultPool.getDataByRef(ctx.data.antigenColumnRef)?.data;
    const values = data?.getDataAsJson<Record<string, string>>()?.["data"];
    if (!values) return undefined;
    return [...new Set(Object.values(values))].sort();
  })

  .output("concentrationUnitLabel", (ctx) => {
    if (ctx.data.concentrationColumnRef === undefined) return undefined;
    const spec = ctx.resultPool.getPColumnSpecByRef(ctx.data.concentrationColumnRef);
    return spec?.annotations?.["pl7.app/label"];
  })

  .output("validationWarnings", (ctx): ValidationIssue[] => {
    const issues: ValidationIssue[] = [];
    const data = ctx.data;

    if (data.concentrationColumnRef !== undefined) {
      const spec = ctx.resultPool.getPColumnSpecByRef(data.concentrationColumnRef);
      const unitLabel = spec?.annotations?.["pl7.app/label"];
      if (unitLabel && /\s/.test(unitLabel)) {
        issues.push({
          severity: "warning",
          message:
            `Concentration column label "${unitLabel}" contains spaces — ` +
            "the full string becomes the K_D,app unit and renders awkwardly. " +
            'Prefer a bare unit string like "nM" or "µM".',
        });
      }
    }

    if (data.antigenColumnRef === undefined && data.targetAntigen) {
      issues.push({
        severity: "warning",
        message:
          "targetAntigen is set but no antigen column is selected — " +
          "the value is ignored. Select an antigen column or clear the targetAntigen field.",
      });
    }

    if (data.antigenColumnRef !== undefined && !data.targetAntigen) {
      issues.push({
        severity: "error",
        message:
          "An antigen column is selected but targetAntigen is empty — " +
          "choose which antigen to analyse before running.",
      });
    }

    if (data.r2ThresholdFailed > data.r2ThresholdGood) {
      issues.push({
        severity: "error",
        message: "Failed R² threshold must be ≤ Good R² threshold.",
      });
    }

    if (data.nMin >= data.nMax) {
      issues.push({
        severity: "error",
        message: "Hill coefficient nMin must be strictly less than nMax.",
      });
    }

    return issues;
  })

  .output("logHandle", (ctx) => ctx.outputs?.resolve("logHandle")?.getLogHandle())

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .outputWithStatus("summaryTable", (ctx) => {
    const pCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPlDataTableV2(ctx, pCols, ctx.data.tableState);
  })

  .outputWithStatus("titrationCurvesPf", (ctx): PFrameHandle | undefined => {
    const signalCols = ctx.outputs?.resolve("signalPf")?.getPColumns();
    if (signalCols === undefined) return undefined;
    const summaryCols = ctx.outputs?.resolve("summaryPf")?.getPColumns() ?? [];
    return createPFrameForGraphs(ctx, [...signalCols, ...summaryCols]);
  })

  .output("titrationCurvesPfCols", (ctx) => {
    const signalCols = ctx.outputs?.resolve("signalPf")?.getPColumns();
    if (signalCols === undefined) return undefined;
    const summaryCols = ctx.outputs?.resolve("summaryPf")?.getPColumns() ?? [];
    return [...signalCols, ...summaryCols].map(
      (c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec,
    );
  })

  .outputWithStatus("summaryPfHandle", (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPFrameForGraphs(ctx, pCols);
  })

  .output("summaryPfCols", (ctx) => {
    const pCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols.map((c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec);
  })

  .output("isEmpty", (ctx) => {
    const pCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols.length === 0;
  })

  .output("binMode", (ctx) => ctx.data.binColumnRef !== undefined)

  .title((ctx) =>
    ctx.data.targetAntigen ? `Titeseq Analysis (${ctx.data.targetAntigen})` : "Titeseq Analysis",
  )

  .subtitle((ctx) => ctx.data.customBlockLabel || ctx.data.defaultBlockLabel)

  .sections((_ctx) => [
    { type: "link", href: "/", label: "Table" },
    { type: "link", href: "/titration-curves", label: "Titration Curves" },
    { type: "link", href: "/kd-distribution", label: "K_D Distribution" },
    { type: "link", href: "/affinity-vs-fit", label: "Affinity vs Fit Quality" },
  ])

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
