import type { GraphMakerState } from "@milaboratories/graph-maker";
import type {
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlDataTableStateV2,
  PlRef,
} from "@platforma-sdk/model";
import {
  BlockModel,
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

export type UiState = {
  tableState: PlDataTableStateV2;
  graphStateTitrationCurves: GraphMakerState;
  graphStateKDHistogram: GraphMakerState;
  graphStateAffinityVsFit: GraphMakerState;
  settingsOpen: boolean;
};

function getDefaultBlockArgs(): BlockArgs {
  return {
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
  };
}

function isNumericValueType(vt: string | undefined): boolean {
  return vt === "Int" || vt === "Long" || vt === "Float" || vt === "Double";
}

function isIntegerValueType(vt: string | undefined): boolean {
  return vt === "Int" || vt === "Long";
}

export const model = BlockModel.create()

  .withArgs<BlockArgs>(getDefaultBlockArgs())

  .withUiState<UiState>({
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
    settingsOpen: true,
  })

  // Block-run gate (R1, R2, R4 blocking checks; numeric validity).
  .argsValid((ctx) => {
    const args = ctx.args;
    if (args.abundanceRef === undefined) return false;
    if (args.concentrationColumnRef === undefined) return false;
    if (args.antigenColumnRef !== undefined && !args.targetAntigen) return false;
    if (args.r2ThresholdFailed > args.r2ThresholdGood) return false;
    if (args.nMin >= args.nMax) return false;
    if (args.r2ThresholdGood < 0 || args.r2ThresholdGood > 1) return false;
    if (args.r2ThresholdFailed < 0 || args.r2ThresholdFailed > 1) return false;
    if (args.minReadsPerConcentration < 1) return false;
    if (args.minConcentrationPoints < 3) return false;
    return true;
  })

  // Abundance PColumn (R1): `[sampleId][clonotypeKey] → numeric`.
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

  // Concentration metadata (R2): per-sample numeric.
  .output("concentrationOptions", (ctx) =>
    ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        isNumericValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    ),
  )

  // Bin metadata (R3): per-sample integer, optional.
  .output("binOptions", (ctx) =>
    ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        isIntegerValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    ),
  )

  // Antigen metadata (R4): per-sample string, optional.
  .output("antigenOptions", (ctx) =>
    ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        spec.valueType === "String" &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    ),
  )

  // Concentration-column unit label (propagated to the `pl7.app/unit`
  // annotation on the `kd` output column in the workflow, R2).
  .output("concentrationUnitLabel", (ctx) => {
    if (ctx.args.concentrationColumnRef === undefined) return undefined;
    const spec = ctx.resultPool.getPColumnSpecByRef(ctx.args.concentrationColumnRef);
    return spec?.annotations?.["pl7.app/label"];
  })

  // Surfaced in the settings drawer as an alert bar.
  .output("validationWarnings", (ctx): ValidationIssue[] => {
    const issues: ValidationIssue[] = [];
    const args = ctx.args;

    if (args.concentrationColumnRef !== undefined) {
      const spec = ctx.resultPool.getPColumnSpecByRef(args.concentrationColumnRef);
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

    if (args.antigenColumnRef === undefined && args.targetAntigen) {
      issues.push({
        severity: "warning",
        message:
          "targetAntigen is set but no antigen column is selected — " +
          "the value is ignored. Select an antigen column or clear the targetAntigen field.",
      });
    }

    if (args.antigenColumnRef !== undefined && !args.targetAntigen) {
      issues.push({
        severity: "error",
        message:
          "An antigen column is selected but targetAntigen is empty — " +
          "choose which antigen to analyse before running.",
      });
    }

    if (args.r2ThresholdFailed > args.r2ThresholdGood) {
      issues.push({
        severity: "error",
        message: "Failed R² threshold must be ≤ Good R² threshold.",
      });
    }

    if (args.nMin >= args.nMax) {
      issues.push({
        severity: "error",
        message: "Hill coefficient nMin must be strictly less than nMax.",
      });
    }

    return issues;
  })

  // Log handle for the Python fitting script (stderr captures WARN lines).
  .output("logHandle", (ctx) => ctx.outputs?.resolve("logHandle")?.getLogHandle())

  // Detects the "workflow running" state for auto-closing the settings drawer.
  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  // Per-clonotype summary columns for the Table tab (R18).
  .outputWithStatus("summaryTable", (ctx) => {
    const pCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPlDataTableV2(ctx, pCols, ctx.uiState.tableState);
  })

  // PFrame for the Titration Curves tab (R15): meanBin scatter +
  // fittedMeanBin curve overlay, faceted by clonotypeKey.
  .outputWithStatus("titrationCurvesPf", (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve("signalPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return createPFrameForGraphs(ctx, pCols);
  })

  // Column specs for default graph options on the Titration Curves page.
  .output("titrationCurvesPfCols", (ctx) => {
    const pCols = ctx.outputs?.resolve("signalPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols.map((c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec);
  })

  // PFrame for the K_D Distribution (R16) and Affinity vs Fit (R17) tabs.
  // Both tabs derive their defaults from the per-clonotype summary columns.
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

  // Empty/failed-state banner signal (R19b, T11). Derived client-side from
  // the summary specs because the workflow currently writes no explicit flag.
  .output("isEmpty", (ctx) => {
    const pCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols.length === 0;
  })

  // "bin" mode is inferred from the presence of a bin column. Used by the
  // UI to swap labels (R19b: "Mean bin" → "Clonotype frequency") and to
  // show the no-bin-mode persistent warning banner.
  .output("binMode", (ctx) => ctx.args.binColumnRef !== undefined)

  .title((ctx) => {
    const parts = ["Titeseq Analysis"];
    const abundanceRef = ctx.args.abundanceRef;
    if (abundanceRef !== undefined) {
      const label =
        ctx.resultPool.getPColumnSpecByRef(abundanceRef)?.annotations?.["pl7.app/label"];
      if (label) parts.push(`— ${label}`);
    }
    if (ctx.args.targetAntigen) parts.push(`(${ctx.args.targetAntigen})`);
    return parts.join(" ");
  })

  .subtitle((ctx) => ctx.args.customBlockLabel || ctx.args.defaultBlockLabel)

  .sections((_ctx) => [
    { type: "link" as const, href: "/" as const, label: "Titration Curves" },
    { type: "link" as const, href: "/kd-distribution" as const, label: "K_D Distribution" },
    { type: "link" as const, href: "/affinity-vs-fit" as const, label: "Affinity vs Fit Quality" },
    { type: "link" as const, href: "/table" as const, label: "Table" },
  ])

  .done(2);

export type BlockOutputs = InferOutputsType<typeof model>;
