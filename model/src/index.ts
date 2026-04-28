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
  plRefsEqual,
} from "@platforma-sdk/model";

export type ValidationIssue = {
  severity: "warning" | "error";
  message: string;
};

export type BlockArgs = {
  abundanceRef?: PlRef;
  concentrationColumnRef?: PlRef;
  binColumnRef?: PlRef;
  antigenColumnRef?: PlRef;
  sortFractionColumnRef?: PlRef;
  targetAntigen?: string;
  minReadsPerConcentration: number;
  minConcentrationPoints: number;
  r2ThresholdGood: number;
  r2ThresholdFailed: number;
  nMin: number;
  nMax: number;
  hookEffectThresholdBin: number;
  hookEffectThresholdNoBin: number;
  hookEffectMinReads: number;
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

// Bind by resolvePath, not block-id domain, so defaults survive a fresh block
// or project. The format — double-encoded JSON — is what GraphMaker stores.
const summaryCol = (name: string, type: "Double" | "String") =>
  `{"kind":"column","name":"{\\"name\\":\\"${name}\\",\\"resolvePath\\":[\\"main\\",\\"summaryPf\\"]}","type":"${type}"}`;

const KD_PLOT_COL = summaryCol("kdPlotPosition", "Double");
const HILL_PLOT_COL = summaryCol("hillPlotPosition", "Double");
const HILL_COEF_COL = summaryCol("hillCoefficient", "Double");
const AFFINITY_CLASS_COL = summaryCol("affinityClass", "String");
const FIT_FAILURE_REASON_COL = summaryCol("fitFailureReason", "String");

// Pin colours by colorIdx so low_r2 stays green and n_out_of_range stays
// purple across datasets. The "light" palette would otherwise assign by
// category-encounter order, which flips when a dataset has only one reason.
const AFFINITY_VS_FIT_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "dots",
  currentTab: null,
  layersSettings: {
    dots: {
      dotFill: { type: "grouping", value: FIT_FAILURE_REASON_COL },
    },
  },
  axesSettings: {
    axisX: { scale: "log" },
  },
  optionsState: {
    type: "scatterplot",
    components: {
      x: { type: "simple", selectorStates: [{ selectedSource: KD_PLOT_COL }] },
      y: { type: "simple", selectorStates: [{ selectedSource: HILL_PLOT_COL }] },
      filters: {
        type: "filter",
        selectorStates: [
          { selectedSource: HILL_COEF_COL, type: "range", selectedFilterRange: { min: 0 } },
          {
            selectedSource: AFFINITY_CLASS_COL,
            type: "equals",
            selectedFilterValues: ["Failed"],
          },
        ],
      },
      grouping: {
        type: "simple",
        selectorStates: [{ selectedSource: FIT_FAILURE_REASON_COL }],
      },
    },
    dividedAxes: {},
  },
  dataBindAes: {
    [FIT_FAILURE_REASON_COL]: {
      type: "categorical",
      palette: "light",
      naAes: { color: "#ccc", lineShape: "solid", dotShape: "21" },
      order: ["low_r2", "n_out_of_range"],
      hidden: {},
      mapping: {
        low_r2: { colorIdx: 0, aes: { color: "#99E099", lineShape: "solid", dotShape: "21" } },
        n_out_of_range: {
          colorIdx: 1,
          aes: { color: "#C1ADFF", lineShape: "solid", dotShape: "21" },
        },
      },
    },
  },
} as GraphMakerState;

const dataModel = new DataModelBuilder().from<BlockData>("v1").init(() => ({
  abundanceRef: undefined,
  concentrationColumnRef: undefined,
  binColumnRef: undefined,
  antigenColumnRef: undefined,
  sortFractionColumnRef: undefined,
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
  defaultBlockLabel: "Tite-Seq Analysis",
  customBlockLabel: "",
  tableState: createPlDataTableStateV2(),
  graphStateTitrationCurves: {
    title: "",
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
    title: "",
    template: "bins",
    currentTab: null,
    layersSettings: {
      bins: { fillColor: "#99E099" },
    },
    axesSettings: {
      axisX: { scale: "log" },
      axisY: { scale: "log" },
      other: { binsCount: 30 },
    },
  },
  graphStateAffinityVsFit: AFFINITY_VS_FIT_DEFAULT_STATE,
  settingsOpen: false,
}));

function isIntegerValueType(vt: string | undefined): boolean {
  return vt === "Int" || vt === "Long";
}

function isFloatValueType(vt: string | undefined): boolean {
  return vt === "Float" || vt === "Double";
}

// Filter per-sample numeric option lists:
//   - drop candidates matching `excludeRef` (cross-exclusion between pickers);
//   - drop all-null columns (they arrive in Python as empty-String and fail validation);
//   - optionally apply `valuePredicate` to the finite values (e.g. [0, 1] range).
// Generic S/O so TS infers the spec type from `ctx.resultPool.getOptions` (normally
// PObjectSpec) without this helper taking a hard dependency on SDK internals.
function filterPopulatedOptions<O extends { ref: PlRef }, S>(
  ctx: {
    resultPool: {
      getOptions: (pred: (spec: S) => boolean) => O[];
      getDataByRef: (
        ref: PlRef,
      ) => { data?: { getDataAsJson: <T>() => T | undefined } } | undefined;
    };
  },
  specPred: (spec: S) => boolean,
  opts: { excludeRef?: PlRef; valuePredicate?: (finite: number[]) => boolean } = {},
): O[] {
  return ctx.resultPool.getOptions(specPred).filter((opt) => {
    if (opts.excludeRef && plRefsEqual(opt.ref, opts.excludeRef)) return false;
    const data = ctx.resultPool.getDataByRef(opt.ref)?.data;
    if (!data) return true;
    const values = data.getDataAsJson<Record<string, number | null>>()?.["data"];
    if (!values) return true;
    const finite = Object.values(values).filter((v): v is number => v !== null && v !== undefined);
    if (finite.length === 0) return false;
    return opts.valuePredicate ? opts.valuePredicate(finite) : true;
  });
}

export const model = BlockModelV3.create(dataModel)

  .args<BlockArgs>((data) => {
    if (data.abundanceRef === undefined) throw new Error("Abundance column is required");
    if (data.concentrationColumnRef === undefined)
      throw new Error("Concentration column is required");
    if (data.antigenColumnRef !== undefined && !data.targetAntigen)
      throw new Error("Target antigen is required when an antigen column is selected");
    if (data.sortFractionColumnRef !== undefined && data.binColumnRef === undefined)
      throw new Error(
        "FACS sort fraction requires a FACS bin column — the correction only applies in bin mode",
      );

    // Reject undefined on the numeric tuning params before any range comparisons
    // run. Without these guards, `undefined < 1` (and friends) silently returns
    // false, so a stored block missing one of these fields would pass validation
    // and produce args with undefined values — breaking the Python fit at runtime.
    if (data.minReadsPerConcentration === undefined)
      throw new Error("Min reads per concentration is required");
    if (data.minConcentrationPoints === undefined)
      throw new Error("Min concentration points is required");
    if (data.r2ThresholdGood === undefined) throw new Error("R² threshold (Good) is required");
    if (data.r2ThresholdFailed === undefined) throw new Error("R² threshold (Failed) is required");
    if (data.nMin === undefined) throw new Error("Hill coefficient nMin is required");
    if (data.nMax === undefined) throw new Error("Hill coefficient nMax is required");
    // Only require the active mode's hook threshold — the other field is hidden
    // in the UI (v-if on binMode), so requiring both would block users from
    // opening Inputs to fix it.
    const isBinMode = data.binColumnRef !== undefined;
    if (isBinMode && data.hookEffectThresholdBin === undefined)
      throw new Error("Hook effect signal-drop threshold (bin mode) is required");
    if (!isBinMode && data.hookEffectThresholdNoBin === undefined)
      throw new Error("Hook effect signal-drop threshold (frequency mode) is required");
    if (data.hookEffectMinReads === undefined)
      throw new Error("Min reads for hook check is required");

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
      sortFractionColumnRef: data.sortFractionColumnRef,
      targetAntigen: data.targetAntigen,
      minReadsPerConcentration: data.minReadsPerConcentration,
      minConcentrationPoints: data.minConcentrationPoints,
      r2ThresholdGood: data.r2ThresholdGood,
      r2ThresholdFailed: data.r2ThresholdFailed,
      nMin: data.nMin,
      nMax: data.nMax,
      // Default the inactive mode's threshold so the args type is satisfied
      // — the workflow only reads the field for the active mode.
      hookEffectThresholdBin: data.hookEffectThresholdBin ?? 0.2,
      hookEffectThresholdNoBin: data.hookEffectThresholdNoBin ?? 0.02,
      hookEffectMinReads: data.hookEffectMinReads,
      // Strings have UI defaults but old stored data may lack them; coerce to
      // empty/safe values rather than throw — they don't break the workflow.
      defaultBlockLabel: data.defaultBlockLabel ?? "Tite-Seq Analysis",
      customBlockLabel: data.customBlockLabel ?? "",
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

  .output("concentrationOptions", (ctx) =>
    filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        isFloatValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
      { excludeRef: ctx.data.sortFractionColumnRef },
    ),
  )

  .output("binOptions", (ctx) =>
    filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        isIntegerValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    ),
  )

  // Disambiguate from concentration by data, not column name: cross-exclude
  // the picked concentration column, then require every value in [0, 1].
  // Molar concentrations (all < 1) still need the cross-exclude to clear.
  .output("sortFractionOptions", (ctx) =>
    filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        isFloatValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
      {
        excludeRef: ctx.data.concentrationColumnRef,
        valuePredicate: (finite) => finite.every((v) => v >= 0 && v <= 1),
      },
    ),
  )

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

    // Surface backend field-level errors (CID conflicts, exec failures, spec
    // validation failures) as severity:"error" ValidationIssues so the user
    // sees them alongside the UI state. Without this, a CID conflict only
    // appears in server logs while the block silently shows empty panels.
    // CIDConflictError is typically transient — the platform's "zebra"
    // self-recovery retries the write — so we label it as such.
    const outputKeys = ["summaryPf", "signalPf", "logHandle"] as const;
    for (const key of outputKeys) {
      const errorAcc = ctx.outputs?.resolve(key)?.getError();
      if (!errorAcc) continue;
      const typeName = errorAcc.resourceType.name;
      const detail = errorAcc.getDataAsJson<{ message?: string }>()?.message ?? typeName;
      const isCidConflict = typeName === "CIDConflictError";
      issues.push({
        severity: "error",
        message: isCidConflict
          ? `Transient result-hash conflict on "${key}": ${detail}. ` +
            `The platform retries automatically; usually self-recovers. ` +
            `If this persists across retries, the workflow likely has a non-deterministic step — report it.`
          : `Workflow output "${key}" failed: ${detail}`,
      });
    }

    if (data.concentrationColumnRef !== undefined) {
      const spec = ctx.resultPool.getPColumnSpecByRef(data.concentrationColumnRef);
      const unitLabel = spec?.annotations?.["pl7.app/label"];
      if (unitLabel && /\s/.test(unitLabel)) {
        issues.push({
          severity: "warning",
          message:
            `Concentration column label "${unitLabel}" contains spaces. ` +
            `The full label becomes the Kd,app unit, so values render as "Kd,app (${unitLabel})". ` +
            'Use a bare unit like "nM" or "µM".',
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

    // Empty targetAntigen surfaces inline: the dropdown's `:error-status` prop
    // shows a red border and `required` shows the asterisk. The args() throw
    // above keeps Run disabled until it's set.

    if (data.sortFractionColumnRef !== undefined && data.binColumnRef === undefined) {
      issues.push({
        severity: "error",
        message:
          "FACS sort fraction requires a FACS bin column — " +
          "the correction only applies in bin mode. Select a FACS bin or clear the sort fraction.",
      });
    }

    // Numeric-field bound violations surface inline on each PlNumberField via
    // its `:error-message` prop (see ui/src/composables/useFieldValidation.ts).
    // The args() throw above gates Run.

    return issues;
  })

  .output("logHandle", (ctx) => ctx.outputs?.resolve("logHandle")?.getLogHandle())

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .outputWithStatus("summaryTable", (ctx) => {
    const summaryCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (summaryCols === undefined) return undefined;
    const signalCols = ctx.outputs?.resolve("signalPf")?.getPColumns() ?? [];

    // Reveal fitFailureReason and the signal columns in this block's Table so
    // users see why each clonotype failed and can export the per-concentration
    // data. Both carry pl7.app/table/visibility: "hidden" to stay out of
    // downstream pickers — overridden locally only.
    const withVisibility =
      (visibility: string) =>
      <T extends { spec: { annotations?: Record<string, string> } }>(c: T): T => ({
        ...c,
        spec: {
          ...c.spec,
          annotations: { ...c.spec.annotations, "pl7.app/table/visibility": visibility },
        },
      });

    const visibleSummary = summaryCols.map((c) =>
      c.spec.name === "pl7.app/vdj/fitFailureReason" ? withVisibility("default")(c) : c,
    );
    const visibleSignal = signalCols.map(withVisibility("default"));

    const kdCol = visibleSummary.find((c) => c.spec.name === "pl7.app/vdj/kd");
    if (!kdCol) return undefined;

    // Include result pool columns whose axes are a direct subset of the anchor
    // (kd) axes — equivalent to enrichment mode / maxHops: 0 in V3.
    // Excludes this block's own outputs (handled above) and File columns.
    const anchorAxisNames = new Set(kdCol.spec.axesSpec.map((a) => a.name));
    const resultPoolCols = (
      ctx.resultPool.getAnchoredPColumns(
        { main: kdCol.spec },
        (spec) =>
          (spec.valueType as string) !== "File" &&
          !spec.annotations?.["pl7.app/trace"]?.includes("milaboratories.titeseq-analysis") &&
          spec.axesSpec.every((a) => anchorAxisNames.has(a.name)),
        { dontWaitAllData: true },
      ) ?? []
    ).map(withVisibility("optional"));

    return createPlDataTableV2(
      ctx,
      [...visibleSummary, ...visibleSignal, ...resultPoolCols],
      ctx.data.tableState,
    );
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

  .title(() => "Tite-Seq Analysis")

  .subtitle((ctx) => ctx.data.customBlockLabel || ctx.data.defaultBlockLabel || "")

  .sections((_ctx) => [
    { type: "link", href: "/", label: "Table" },
    { type: "link", href: "/titration-curves", label: "Titration Curves" },
    { type: "link", href: "/kd-distribution", label: "Kd Distribution" },
    { type: "link", href: "/affinity-vs-fit", label: "Affinity vs Fit Quality" },
  ])

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
