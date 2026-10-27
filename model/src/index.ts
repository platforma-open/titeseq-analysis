import type { GraphMakerState } from "@milaboratories/graph-maker";
import type {
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlDataTableStateV2,
  PlRef,
} from "@platforma-sdk/model";
import type { BlockParams, FitModel } from "@platforma-open/platforma-open.titeseq-analysis.kind";
import { kind } from "@platforma-open/platforma-open.titeseq-analysis.kind";
import { deriveTemplateParams } from "./templateParams";

export type { BlockParams, FitModel } from "@platforma-open/platforma-open.titeseq-analysis.kind";
export { deriveTemplateParams };
import {
  BlockModelV3,
  ColumnsCollection,
  DataColumn,
  DataModelBuilder,
  createPFrameForGraphs,
  createPlDataTableStateV2,
  createPlDataTableV3,
  isPColumnSpec,
  plRefsEqual,
} from "@platforma-sdk/model";

export type ValidationIssue = {
  severity: "warning" | "error";
  message: string;
};

/**
 * What the workflow consumes. `fitModel` here is the model the block RESOLVED, not a stored
 * choice: `deriveArgs` routes it from the bound columns unless `fitModelOverride` forces one.
 */
export type BlockArgs = {
  abundanceRef?: PlRef;
  concentrationColumnRef?: PlRef;
  binColumnRef?: PlRef;
  antigenColumnRef?: PlRef;
  sortFractionColumnRef?: PlRef;
  // Cells sorted into each sample (C_bc in Adams 2016 eq. E4). Preferred over sort fractions: it
  // also lets the fit cap a bin's reads at its cell count.
  cellCountColumnRef?: PlRef;
  // Measured mean log10 fluorescence of the cells in each sample's gate. Bound, each bin's
  // fluorescence is read per concentration instead of from one fixed gate spacing.
  binFluorescenceColumnRef?: PlRef;
  // FACS gate width in decades, per sample. Bound, it supplies the bin width.
  gateWidthColumnRef?: PlRef;
  targetAntigen?: string;
  fitModel: FitModel;
  binWidthDecades: number;
  minPlateauFraction: number;
  baselineFluorescence?: number;
  amplitudePriorQ?: number;
  hookDetectorEnabled: boolean;
  pinHillExponent: boolean;
  kdSdGood: number;
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

/** What the settings drawer edits and the project stores. */
export type BlockData = Omit<BlockArgs, "fitModel"> & {
  // Forces the curve model instead of routing on the bound columns. Undefined means route.
  fitModelOverride?: FitModel;
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

// Pin colours by colorIdx so convergence_failure stays green and
// insufficient_dynamic_range stays purple across datasets. The "light" palette would otherwise
// assign by category-encounter order, which flips when a dataset has only one reason.
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
      order: ["convergence_failure", "insufficient_dynamic_range"],
      hidden: {},
      mapping: {
        convergence_failure: {
          colorIdx: 0,
          aes: { color: "#99E099", lineShape: "solid", dotShape: "21" },
        },
        insufficient_dynamic_range: {
          colorIdx: 1,
          aes: { color: "#C1ADFF", lineShape: "solid", dotShape: "21" },
        },
      },
    },
  },
} as GraphMakerState;

// Adams 2016 used gate widths of 0.61-0.71 decades and plateau fractions of 0.63-0.70.
const DEFAULT_BIN_WIDTH_DECADES = 0.65;
const DEFAULT_MIN_PLATEAU_FRACTION = 0.65;
const DEFAULT_KD_SD_GOOD = 0.5;
const DEFAULT_HOOK_DETECTOR_ENABLED = false;
const DEFAULT_PIN_HILL_EXPONENT = true;

/** v2: the Adams 2016 model, chosen by hand, with the uncorrected-fit escape hatch. */
export type BlockDataV2 = Omit<
  BlockData,
  | "fitModelOverride"
  | "binFluorescenceColumnRef"
  | "gateWidthColumnRef"
  | "amplitudePriorQ"
  | "hookDetectorEnabled"
  | "pinHillExponent"
> & {
  fitModel: FitModel;
  allowUncorrected: boolean;
};

/** Block data before the Adams 2016 model: no cell counts, no model choice. */
export type BlockDataV1 = Omit<
  BlockDataV2,
  | "cellCountColumnRef"
  | "fitModel"
  | "binWidthDecades"
  | "minPlateauFraction"
  | "baselineFluorescence"
  | "allowUncorrected"
  | "kdSdGood"
>;

/**
 * A fresh block's data, seeded by whatever the creator or template supplied. Every param the
 * contract carries is honoured here, and every field it does not carry falls back to the
 * block's own default — so a block created without params is exactly the block this returns.
 */
export const initBlockData = (params?: BlockParams): BlockData => ({
  abundanceRef: params?.abundanceRef,
  concentrationColumnRef: params?.concentrationColumnRef,
  binColumnRef: params?.binColumnRef,
  antigenColumnRef: params?.antigenColumnRef,
  sortFractionColumnRef: params?.sortFractionColumnRef,
  cellCountColumnRef: params?.cellCountColumnRef,
  binFluorescenceColumnRef: params?.binFluorescenceColumnRef,
  gateWidthColumnRef: params?.gateWidthColumnRef,
  targetAntigen: params?.targetAntigen,
  fitModelOverride: params?.fitModelOverride,
  binWidthDecades: params?.binWidthDecades ?? DEFAULT_BIN_WIDTH_DECADES,
  minPlateauFraction: params?.minPlateauFraction ?? DEFAULT_MIN_PLATEAU_FRACTION,
  baselineFluorescence: params?.baselineFluorescence,
  amplitudePriorQ: params?.amplitudePriorQ,
  hookDetectorEnabled: params?.hookDetectorEnabled ?? DEFAULT_HOOK_DETECTOR_ENABLED,
  pinHillExponent: params?.pinHillExponent ?? DEFAULT_PIN_HILL_EXPONENT,
  kdSdGood: params?.kdSdGood ?? DEFAULT_KD_SD_GOOD,
  minReadsPerConcentration: params?.minReadsPerConcentration ?? 3,
  minConcentrationPoints: params?.minConcentrationPoints ?? 5,
  r2ThresholdGood: params?.r2ThresholdGood ?? 0.8,
  r2ThresholdFailed: params?.r2ThresholdFailed ?? 0.5,
  nMin: params?.nMin ?? 0.5,
  nMax: params?.nMax ?? 2.0,
  hookEffectThresholdBin: params?.hookEffectThresholdBin ?? 0.2,
  hookEffectThresholdNoBin: params?.hookEffectThresholdNoBin ?? 0.02,
  hookEffectMinReads: params?.hookEffectMinReads ?? 20,
  defaultBlockLabel: "Tite-Seq Analysis",
  customBlockLabel: params?.customBlockLabel ?? "",
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
});

/**
 * v1 -> v2. A project configured before v2 has no cell counts; it stays on the linear Hill it ran
 * with, so a rerun reproduces its results instead of failing.
 */
export const migrateV1ToV2 = (v1: BlockDataV1): BlockDataV2 => ({
  ...v1,
  cellCountColumnRef: undefined,
  fitModel: "hill_free_n",
  binWidthDecades: DEFAULT_BIN_WIDTH_DECADES,
  minPlateauFraction: DEFAULT_MIN_PLATEAU_FRACTION,
  baselineFluorescence: undefined,
  allowUncorrected: false,
  kdSdGood: DEFAULT_KD_SD_GOOD,
});

/**
 * v2 -> v3. Drops the stored model choice so every project is routed on its bound columns, and
 * drops the uncorrected-fit flag, which no longer gates anything.
 *
 * The pin has to be dropped rather than carried: `migrate` fast-paths a project already at its
 * latest version, so a project that passed through v1 -> v2 would otherwise keep
 * `fitModel: "hill_free_n"` forever and the routing would never reach it. Clearing it here also
 * keeps a legacy pin distinguishable from a deliberate override, which shares the field otherwise.
 */
export const migrateV2ToV3 = ({
  fitModel: _fitModel,
  allowUncorrected: _allowUncorrected,
  ...v2
}: BlockDataV2): BlockData => ({
  ...v2,
  fitModelOverride: undefined,
  binFluorescenceColumnRef: undefined,
  gateWidthColumnRef: undefined,
  amplitudePriorQ: undefined,
  hookDetectorEnabled: DEFAULT_HOOK_DETECTOR_ENABLED,
  pinHillExponent: DEFAULT_PIN_HILL_EXPONENT,
});

const dataModel = new DataModelBuilder({ kind })
  .from<BlockDataV1>("v1")
  .migrate<BlockDataV2>("v2", migrateV1ToV2)
  .migrate<BlockData>("v3", migrateV2ToV3)
  .init(({ params }) => initBlockData(params));

function isIntegerValueType(vt: string | undefined): boolean {
  return vt === "Int" || vt === "Long";
}

// The Python cause of a failed run, from the error chain an errored output throws.
function runFailureCause(e: unknown): string {
  const text = e instanceof Error ? e.message : String(e);
  const cause = text.match(/InputValidationError: ([^\n]+)/);
  if (cause) return cause[1].trim();
  const lines = text.split("\n").filter((l) => l.trim() !== "");
  return (lines[lines.length - 1] ?? "unknown error").trim().slice(0, 500);
}

function isFloatValueType(vt: string | undefined): boolean {
  return vt === "Float" || vt === "Double";
}

// Column selectors read a bare string as a regex, so literal names and domain
// values go through an explicit exact matcher.
function exactly(value: string) {
  return { type: "exact" as const, value };
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

type HillSettings = Pick<
  BlockArgs,
  | "r2ThresholdGood"
  | "r2ThresholdFailed"
  | "nMin"
  | "nMax"
  | "hookEffectThresholdBin"
  | "hookEffectThresholdNoBin"
  | "hookEffectMinReads"
>;

const hillSettingsOf = (data: HillSettings): HillSettings => ({
  r2ThresholdGood: data.r2ThresholdGood,
  r2ThresholdFailed: data.r2ThresholdFailed,
  nMin: data.nMin,
  nMax: data.nMax,
  hookEffectThresholdBin: data.hookEffectThresholdBin,
  hookEffectThresholdNoBin: data.hookEffectThresholdNoBin,
  hookEffectMinReads: data.hookEffectMinReads,
});

const HILL_SETTINGS_UNUSED: HillSettings = hillSettingsOf(initBlockData());

const validHillSettings = (data: BlockData, isBinMode: boolean): HillSettings => {
  if (data.r2ThresholdGood === undefined) throw new Error("R² threshold (Good) is required");
  if (data.r2ThresholdFailed === undefined) throw new Error("R² threshold (Failed) is required");
  if (data.nMin === undefined) throw new Error("Hill coefficient nMin is required");
  if (data.nMax === undefined) throw new Error("Hill coefficient nMax is required");
  // Only require the active mode's hook threshold — the other field is hidden
  // in the UI (v-if on binMode), so requiring both would block users from
  // opening Inputs to fix it.
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
  return {
    ...hillSettingsOf(data),
    // Default the inactive mode's threshold so the args type is satisfied
    // — the workflow only reads the field for the active mode.
    hookEffectThresholdBin: data.hookEffectThresholdBin ?? 0.2,
    hookEffectThresholdNoBin: data.hookEffectThresholdNoBin ?? 0.02,
  };
};

/**
 * The curve model the bound columns support.
 *
 * Measurement needs the cells sorted into each sample. Without them the mean bin is uncorrected,
 * and feeding that to the measurement likelihood is measurably WORSE than the ranking fit
 * (rank agreement 0.888 -> 0.822 on Adams 2016 replicate 1). "More data, better results" is a
 * property of this routing, not of the inputs: more data has to switch the model to hold.
 */
export const routeFitModel = (
  data: Pick<BlockData, "binColumnRef" | "cellCountColumnRef" | "sortFractionColumnRef">,
): FitModel => {
  if (data.binColumnRef === undefined) return "hill_free_n";
  const cellsKnown =
    data.cellCountColumnRef !== undefined || data.sortFractionColumnRef !== undefined;
  return cellsKnown ? "poisson_log" : "hill_free_n";
};

export const deriveArgs = (data: BlockData): BlockArgs => {
  if (data.abundanceRef === undefined) throw new Error("Abundance column is required");
  if (data.concentrationColumnRef === undefined)
    throw new Error("Concentration column is required");
  if (data.antigenColumnRef !== undefined && !data.targetAntigen)
    throw new Error("Target antigen is required when an antigen column is selected");
  if (data.sortFractionColumnRef !== undefined && data.binColumnRef === undefined)
    throw new Error(
      "FACS sort fraction requires a FACS bin column — the correction only applies in bin mode",
    );
  if (data.cellCountColumnRef !== undefined && data.binColumnRef === undefined)
    throw new Error("Sorted cell counts require a FACS bin column — they only apply in bin mode");
  if (data.binFluorescenceColumnRef !== undefined && data.binColumnRef === undefined)
    throw new Error("Measured bin fluorescence requires a FACS bin column");
  if (data.gateWidthColumnRef !== undefined && data.binColumnRef === undefined)
    throw new Error("FACS gate width requires a FACS bin column");
  if (!(data.binWidthDecades > 0)) throw new Error("Bin width must be above 0 decades");
  if (data.amplitudePriorQ !== undefined && !(data.amplitudePriorQ > 0 && data.amplitudePriorQ < 1))
    throw new Error("Amplitude prior must be between 0 and 1, exclusive");
  if (!(data.minPlateauFraction > 0 && data.minPlateauFraction < 1))
    throw new Error("Minimum plateau must be between 0 and 1, exclusive");
  if (!(data.kdSdGood > 0)) throw new Error("Good-class Kd,app uncertainty must be above 0");
  if (data.baselineFluorescence !== undefined && !(data.baselineFluorescence > 0))
    throw new Error("0 M baseline fluorescence must be above 0");

  // Reject undefined on the numeric tuning params before any range comparisons
  // run. Without these guards, `undefined < 1` (and friends) silently returns
  // false, so a stored block missing one of these fields would pass validation
  // and produce args with undefined values — breaking the Python fit at runtime.
  if (data.minReadsPerConcentration === undefined)
    throw new Error("Min reads per concentration is required");
  if (data.minConcentrationPoints === undefined)
    throw new Error("Min concentration points is required");
  if (data.minReadsPerConcentration < 1) throw new Error("Min reads per concentration must be ≥ 1");
  if (data.minConcentrationPoints < 3) throw new Error("Min concentration points must be ≥ 3");

  // Only the linear Hill fit reads the R², Hill n and hook settings, and the UI hides them under
  // the Adams model. They are validated only when used and sent as fixed values otherwise, so a
  // hidden setting can neither block a run nor mark the block stale.
  const isBinMode = data.binColumnRef !== undefined;
  const fitModel = data.fitModelOverride ?? routeFitModel(data);
  const usesHill = fitModel === "hill_free_n";
  const hill = usesHill ? validHillSettings(data, isBinMode) : HILL_SETTINGS_UNUSED;
  return {
    abundanceRef: data.abundanceRef,
    concentrationColumnRef: data.concentrationColumnRef,
    binColumnRef: data.binColumnRef,
    antigenColumnRef: data.antigenColumnRef,
    sortFractionColumnRef: data.sortFractionColumnRef,
    cellCountColumnRef: data.cellCountColumnRef,
    binFluorescenceColumnRef: data.binFluorescenceColumnRef,
    gateWidthColumnRef: data.gateWidthColumnRef,
    targetAntigen: data.targetAntigen,
    fitModel,
    binWidthDecades: data.binWidthDecades,
    minPlateauFraction: data.minPlateauFraction,
    baselineFluorescence: data.baselineFluorescence,
    amplitudePriorQ: data.amplitudePriorQ,
    hookDetectorEnabled: data.hookDetectorEnabled ?? DEFAULT_HOOK_DETECTOR_ENABLED,
    pinHillExponent: data.pinHillExponent ?? DEFAULT_PIN_HILL_EXPONENT,
    kdSdGood: data.kdSdGood,
    minReadsPerConcentration: data.minReadsPerConcentration,
    minConcentrationPoints: data.minConcentrationPoints,
    ...hill,
    // Strings have UI defaults but old stored data may lack them; coerce to
    // empty/safe values rather than throw — they don't break the workflow.
    defaultBlockLabel: data.defaultBlockLabel ?? "Tite-Seq Analysis",
    customBlockLabel: data.customBlockLabel ?? "",
  };
};

export const platforma = BlockModelV3.create({ dataModel, kind })

  .templateParams(deriveTemplateParams)

  .args<BlockArgs>(deriveArgs)

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
        {
          axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/variantKey" }],
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

  // Cells sorted per sample: any non-negative per-sample number other than the columns already
  // picked as concentration, bin or sort fraction.
  .output("cellCountOptions", (ctx) =>
    filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        (isFloatValueType(spec.valueType) || isIntegerValueType(spec.valueType)) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
      { valuePredicate: (finite) => finite.every((v) => v >= 0) },
    ).filter(
      (opt) =>
        ![
          ctx.data.concentrationColumnRef,
          ctx.data.binColumnRef,
          ctx.data.sortFractionColumnRef,
        ].some((ref) => ref !== undefined && plRefsEqual(opt.ref, ref)),
    ),
  )

  // Measured mean log10 fluorescence of each sample's gate. Any finite per-sample number other
  // than the columns already picked; log fluorescence may be negative, so no value predicate.
  .output("binFluorescenceOptions", (ctx) =>
    filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        (isFloatValueType(spec.valueType) || isIntegerValueType(spec.valueType)) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    ).filter(
      (opt) =>
        ![
          ctx.data.concentrationColumnRef,
          ctx.data.binColumnRef,
          ctx.data.sortFractionColumnRef,
          ctx.data.cellCountColumnRef,
          ctx.data.gateWidthColumnRef,
        ].some((ref) => ref !== undefined && plRefsEqual(opt.ref, ref)),
    ),
  )

  // FACS gate width in decades. Strictly positive and small — a gate spanning more than a decade
  // would make the bin scale meaningless.
  .output("gateWidthOptions", (ctx) =>
    filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        isFloatValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
      { valuePredicate: (finite) => finite.every((v) => v > 0 && v <= 2) },
    ).filter(
      (opt) =>
        ![
          ctx.data.concentrationColumnRef,
          ctx.data.binColumnRef,
          ctx.data.sortFractionColumnRef,
          ctx.data.cellCountColumnRef,
          ctx.data.binFluorescenceColumnRef,
        ].some((ref) => ref !== undefined && plRefsEqual(opt.ref, ref)),
    ),
  )

  // The model the bound columns route to, or the override. Shown beside the results, because a
  // setting chosen for the user rather than by the user has to be visible where the result is read.
  .output(
    "resolvedFitModel",
    (ctx): FitModel => ctx.data.fitModelOverride ?? routeFitModel(ctx.data),
  )

  .output("fitModelIsOverridden", (ctx) => ctx.data.fitModelOverride !== undefined)

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
    // Resolving a field that holds an error and no value throws; the throw carries the cause.
    for (const key of outputKeys) {
      let errorAcc;
      try {
        errorAcc = ctx.outputs?.resolve(key)?.getError();
      } catch (e) {
        issues.push({ severity: "error", message: `The run failed: ${runFailureCause(e)}` });
        break;
      }
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

    // Nothing below blocks a run. Each says what the block did with the data it was given and
    // what one more input would buy.
    const routed = data.fitModelOverride ?? routeFitModel(data);

    if (data.binColumnRef !== undefined && routed === "hill_free_n" && !data.fitModelOverride) {
      issues.push({
        severity: "warning",
        message:
          "No sorted cell counts: the ranking fit ran, which orders candidates but does not " +
          "measure them. Supplying the cells sorted into each sample switches the block to the " +
          "Adams 2016 measurement model and makes Kd,app comparable between runs.",
      });
    }

    if (data.baselineFluorescence === undefined && data.concentrationColumnRef !== undefined) {
      const conc = ctx.resultPool
        .getDataByRef(data.concentrationColumnRef)
        ?.data?.getDataAsJson<Record<string, Record<string, number | null>>>()?.["data"];
      if (conc && !Object.values(conc).some((v) => v === 0)) {
        issues.push({
          severity: "warning",
          message:
            "No 0 M (no-antigen) sample and no baseline entered: the background fluorescence " +
            "floated as a free parameter, so absolute Kd,app values are not comparable between " +
            "runs. Ranking within this run is unaffected. Add a 0 M control, or enter the 0 M " +
            "baseline fluorescence under Curve model.",
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

    const kdCol = summaryCols.find((c) => c.spec.name === "pl7.app/vdj/kd");
    if (!kdCol) return undefined;

    // Include result pool columns keyed by the same clonotype axis as the kd
    // column — i.e., the same cohort of clonotypes. Each axis is matched by
    // full spec (name + domain), not just name, so when the upstream is a
    // redefine-clonotypes run, columns keyed by the original (pre-redefine)
    // clonotype axis are correctly excluded — those have the same axis name but
    // a different domain and a disjoint key space, which would otherwise fail
    // spec integration with "axes sets are disjoint".
    const resultPoolCols = ColumnsCollection(["result_pool"], { ctx: ctx.ctx })
      .filter({
        include: {
          axes: kdCol.spec.axesSpec.map((axis) => ({
            name: [exactly(axis.name)],
            domain: Object.fromEntries(
              Object.entries(axis.domain ?? {}).map(([key, value]) => [key, [exactly(value)]]),
            ),
          })),
          partialAxesMatch: false,
        },
        // The table attaches axis labels itself; a copy re-annotated below
        // would get its own id and show up as a second label column.
        exclude: { name: [exactly("pl7.app/label")] },
      })
      .getColumns()
      .filter((c) => {
        const spec = c.getSpec();
        return (
          (spec.valueType as string) !== "File" &&
          !spec.annotations?.["pl7.app/trace"]?.includes("milaboratories.titeseq-analysis")
        );
      })
      .map((c) => c.withSpecs({ annotations: { "pl7.app/table/visibility": "optional" } }));

    return createPlDataTableV3(ctx, {
      primaryColumns: [...summaryCols, ...signalCols].map((c) => DataColumn.fromColumn(c)),
      columns: resultPoolCols,
      tableState: ctx.data.tableState,
      displayOptions: {
        // Reveal fitFailureReason and the signal columns in this block's Table
        // so users see why each clonotype failed and can export the
        // per-concentration data. Both carry pl7.app/table/visibility: "hidden"
        // to stay out of downstream pickers.
        visibility: [
          {
            match: {
              name: [
                exactly("pl7.app/vdj/fitFailureReason"),
                ...signalCols.map((c) => exactly(c.spec.name)),
              ],
            },
            visibility: "default",
          },
        ],
      },
    });
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

  // True when the run produced no clonotype that is not Failed; undefined before results exist.
  .output("isEmpty", (ctx) => {
    const summary = ctx.outputs
      ?.resolve("fitSummary")
      ?.getDataAsJson<{ clonotypes: number; failed: number }>();
    if (summary === undefined) return undefined;
    return summary.failed === summary.clonotypes;
  })

  // The mode of the results on screen: the args of the last successful run, not unsaved settings.
  .output("binMode", (ctx) =>
    ctx.activeArgs === undefined ? undefined : ctx.activeArgs.binColumnRef !== undefined,
  )

  .title(() => "Tite-Seq Analysis")

  .subtitle((ctx) => ctx.data.customBlockLabel || ctx.data.defaultBlockLabel || "")

  .sections((_ctx) => [
    { type: "link", href: "/", label: "Table" },
    { type: "link", href: "/titration-curves", label: "Titration Curves" },
    { type: "link", href: "/kd-distribution", label: "Kd Distribution" },
    { type: "link", href: "/affinity-vs-fit", label: "Affinity vs Fit Quality" },
  ])

  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
