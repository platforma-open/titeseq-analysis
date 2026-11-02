import type { GraphMakerState } from "@milaboratories/graph-maker";
import type {
  DatasetOption,
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlDataTableStateV2,
  PlRef,
  RenderCtx,
} from "@platforma-sdk/model";
import type {
  BlockParams,
  ReplicateErrorModel,
} from "@platforma-open/platforma-open.titeseq-analysis.kind";
import {
  DEFAULT_CENSORED_MAJORITY,
  DEFAULT_REPLICATE_ERROR_MODEL,
  kind,
} from "@platforma-open/platforma-open.titeseq-analysis.kind";
import type { FitSettings, FitSettingsDraft, PresetId } from "./presets";
import type { LogFluorescenceSetting } from "./presets";
import {
  DEFAULT_EQUATION,
  DEFAULT_PLATEAU_FLOOR,
  DEFAULT_R2_GATE,
  HILL_GATE_INDEX_DEFAULTS,
} from "./presets";
import { deriveTemplateParams } from "./templateParams";

export type {
  BlockParams,
  Equation,
  ReplicateErrorModel,
} from "@platforma-open/platforma-open.titeseq-analysis.kind";
export {
  DEFAULT_CENSORED_MAJORITY,
  DEFAULT_REPLICATE_ERROR_MODEL,
  EQUATIONS,
  REPLICATE_ERROR_MODELS,
} from "@platforma-open/platforma-open.titeseq-analysis.kind";
export * from "./presets";
export { deriveTemplateParams };
import {
  BlockModelV3,
  ColumnsCollection,
  DataColumn,
  DataModelBuilder,
  buildDatasetOptions,
  createGlobalPObjectId,
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

/** The columns a titration is read from, and the antigen it is read for. */
type FitInputs = {
  abundanceRef?: PlRef;
  concentrationColumnRef?: PlRef;
  binColumnRef?: PlRef;
  antigenColumnRef?: PlRef;
  sortFractionColumnRef?: PlRef;
  /**
   * The gate's measured mean log10 fluorescence, one value a (concentration, gate). Read by
   * `log_fluorescence` only, which otherwise assumes the gates sit an equal distance apart.
   */
  binFluorescenceColumnRef?: PlRef;
  targetAntigen?: string;
  /**
   * Which replicate each sample belongs to. Each replicate's titration is fitted on its own and
   * the per-replicate values are then combined. Unset, the block fits one titration.
   */
  replicateColumnRef?: PlRef;
};

/** Every replicate setting the workflow reads. Projected only when a replicate column is set. */
type ReplicateSettings = {
  replicateErrorModel: ReplicateErrorModel;
  /** Decades. Absent leaves the fit-error filter off. */
  replicateFitErrorCutoff?: number;
  /** Absent leaves the R² filter off. */
  replicateR2Floor?: number;
  replicateCensoredMajority: number;
};

type BlockLabels = {
  defaultBlockLabel: string;
  customBlockLabel: string;
};

/** What every equation's args carry: the inputs except gate fluorescence, replicates, labels and
 *  subset. */
type ArgsInputs = Omit<FitInputs, "binFluorescenceColumnRef" | "replicateColumnRef">;
/**
 * Appended last to every equation's args. Every replicate setting is present only beside a
 * replicate column, so a block that set none projects the args it always projected.
 */
type ArgsLabels = BlockLabels & {
  inputFilter?: string;
  replicateColumnRef?: PlRef;
} & Partial<ReplicateSettings>;
type CommonArgs = ArgsInputs & ArgsLabels;

/**
 * Each equation's args carry only the settings that equation reads, so a setting the drawer hides
 * can neither block Run nor stale the block.
 *
 * `hill_gate_index` omits `equation`. The workflow and the Python `FitParams` read an absent
 * `equation` as `hill_gate_index`, so a project created on v2.4.0 or earlier projects the args it
 * always projected, and the upgrade does not stale it.
 */
export type BlockArgs =
  | (CommonArgs & Omit<FitSettings, "equation" | LogFluorescenceSetting> & { equation?: undefined })
  | (CommonArgs & {
      equation: "log_fluorescence";
      plateauFloor: number;
      /** Present only while the gate is on, so an ungated block's args name no gate. */
      r2GateMin?: number;
      /** Present only while the gate is on and an SE limit is set. */
      kdSeGateMax?: number;
      /** Present only as true. */
      requireAllConcentrations?: true;
      binFluorescenceColumnRef?: PlRef;
    });

export type BlockData = FitInputs &
  BlockLabels &
  FitSettingsDraft &
  Pick<
    ReplicateSettings,
    "replicateErrorModel" | "replicateFitErrorCutoff" | "replicateR2Floor"
  > & {
    /** The censoring majority as the drawer edits it: a cleared field reads `undefined`. */
    replicateCensoredMajority?: number;
    /** The per-replicate table, shown only when a run combined replicates. */
    replicateTableState: PlDataTableStateV2;
    graphStateReplicates: GraphMakerState;
    /** Subset column picked with the abundance in PlDatasetSelector (e.g. a Repertoire Labeling label). */
    filterRef?: PlRef;
    /** The preset last applied, kept so the drawer can name a preset the scientist then edited.
     *  Not authoritative: the settings themselves are. Never reaches the workflow. */
    appliedPreset?: PresetId;
    tableState: PlDataTableStateV2;
    graphStateTitrationCurves: GraphMakerState;
    graphStateKDHistogram: GraphMakerState;
    graphStateAffinityVsFit: GraphMakerState;
    settingsOpen: boolean;
  };

/**
 * v1: the block as published through 2.4.0 — one fit, the free-exponent Hill on the read-weighted
 * mean gate index, with no equation to choose. This is the shape every released project stores.
 */
export type BlockDataV1 = Omit<BlockDataV2, "equation" | "plateauFloor" | "appliedPreset">;

/** v2: the equation choice, before the log-fluorescence R² gate. Never released. */
export type BlockDataV2 = Omit<BlockDataV3, "r2GateEnabled" | "r2GateMin">;

/** v3: the R² gate, before sample replicates. Never released. */
export type BlockDataV3 = Omit<
  BlockData,
  | "replicateColumnRef"
  | "replicateErrorModel"
  | "replicateFitErrorCutoff"
  | "replicateR2Floor"
  | "replicateCensoredMajority"
  | "replicateTableState"
  | "graphStateReplicates"
>;

/**
 * The per-replicate chart: one dot per replicate on a log Kd,app axis. The page binds the sources
 * from the columns the run produced.
 */
const REPLICATES_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "dots",
  currentTab: null,
  layersSettings: { dots: {} },
  axesSettings: { axisY: { scale: "log" } },
} as GraphMakerState;

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
    // New blocks use log_fluorescence, whose y is the SE. useFitSettings sets y per equation.
    axisY: { scale: "log" },
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

/**
 * A fresh block's data, seeded by whatever the creator or template supplied. Every param the
 * contract carries is honoured here, and every field it does not carry falls back to the
 * block's own default — so a block created without params is exactly the block this returns.
 */
export const initBlockData = (params?: BlockParams): BlockData => ({
  abundanceRef: params?.abundanceRef,
  filterRef: params?.filterRef,
  concentrationColumnRef: params?.concentrationColumnRef,
  binColumnRef: params?.binColumnRef,
  antigenColumnRef: params?.antigenColumnRef,
  sortFractionColumnRef: params?.sortFractionColumnRef,
  targetAntigen: params?.targetAntigen,
  // A new block fits the log-fluorescence equation. An args projection that omits `equation`
  // means the other one — see `BlockArgs`. So does a template exported before the
  // equation existed: it carries the fit settings but no equation.
  equation:
    params?.equation ??
    (params?.minReadsPerConcentration !== undefined ? "hill_gate_index" : DEFAULT_EQUATION),
  plateauFloor: params?.plateauFloor ?? DEFAULT_PLATEAU_FLOOR,
  r2GateEnabled: params?.r2GateEnabled ?? DEFAULT_R2_GATE.r2GateEnabled,
  r2GateMin: params?.r2GateMin ?? DEFAULT_R2_GATE.r2GateMin,
  replicateColumnRef: params?.replicateColumnRef,
  replicateErrorModel: params?.replicateErrorModel ?? DEFAULT_REPLICATE_ERROR_MODEL,
  // Both filters start off: no rule that is not verified against a paper's own code runs by
  // default. The censoring rule has no paper at all, so its threshold is the block's own.
  replicateFitErrorCutoff: params?.replicateFitErrorCutoff,
  replicateR2Floor: params?.replicateR2Floor,
  replicateCensoredMajority: params?.replicateCensoredMajority ?? DEFAULT_CENSORED_MAJORITY,
  replicateTableState: createPlDataTableStateV2(),
  graphStateReplicates: REPLICATES_DEFAULT_STATE,
  minReadsPerConcentration:
    params?.minReadsPerConcentration ?? HILL_GATE_INDEX_DEFAULTS.minReadsPerConcentration,
  minConcentrationPoints:
    params?.minConcentrationPoints ?? HILL_GATE_INDEX_DEFAULTS.minConcentrationPoints,
  r2ThresholdGood: params?.r2ThresholdGood ?? HILL_GATE_INDEX_DEFAULTS.r2ThresholdGood,
  r2ThresholdFailed: params?.r2ThresholdFailed ?? HILL_GATE_INDEX_DEFAULTS.r2ThresholdFailed,
  nMin: params?.nMin ?? HILL_GATE_INDEX_DEFAULTS.nMin,
  nMax: params?.nMax ?? HILL_GATE_INDEX_DEFAULTS.nMax,
  hookEffectThresholdBin:
    params?.hookEffectThresholdBin ?? HILL_GATE_INDEX_DEFAULTS.hookEffectThresholdBin,
  hookEffectThresholdNoBin:
    params?.hookEffectThresholdNoBin ?? HILL_GATE_INDEX_DEFAULTS.hookEffectThresholdNoBin,
  hookEffectMinReads: params?.hookEffectMinReads ?? HILL_GATE_INDEX_DEFAULTS.hookEffectMinReads,
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
 * v1 -> v2. A v1 project (v2.4.0 or earlier) ran the gate-index fit, so name it.
 *
 * Writing the field is the point. An absent `equation` reads as a new block's default,
 * `log_fluorescence`, and every released project would change the equation it fits.
 *
 * The plateau floor arrives at a new block's value and stays unread. `hill_gate_index` does
 * not read one. Nothing else moves: every setting, reference, label and piece of view state
 * crosses untouched, so a re-run reproduces the numbers the project already holds.
 */
export const migrateV1ToV2 = (v1: BlockDataV1): BlockDataV2 => ({
  ...v1,
  equation: "hill_gate_index",
  plateauFloor: DEFAULT_PLATEAU_FLOOR,
});

/**
 * v2 -> v3. The R² gate arrives off, so a log-fluorescence project rates its fits as it did. v2
 * shipped only in development builds; the step repairs their stored projects.
 */
export const migrateV2ToV3 = (v2: BlockDataV2): BlockDataV3 => ({ ...v2, ...DEFAULT_R2_GATE });

/**
 * v3 -> v4. Sample replicates arrive unset: no replicate column, both filters off, the default
 * error model and censoring majority, and empty per-replicate view state. The args a v3 project
 * sends do not move, so it does not re-run.
 */
export const migrateV3ToV4 = (v3: BlockDataV3): BlockData => ({
  ...v3,
  replicateErrorModel: DEFAULT_REPLICATE_ERROR_MODEL,
  replicateCensoredMajority: DEFAULT_CENSORED_MAJORITY,
  replicateTableState: createPlDataTableStateV2(),
  graphStateReplicates: REPLICATES_DEFAULT_STATE,
});

export const dataModel = new DataModelBuilder({ kind })
  .from<BlockDataV1>("v1")
  .migrate<BlockDataV2>("v2", migrateV1ToV2)
  .migrate<BlockDataV3>("v3", migrateV2ToV3)
  .migrate<BlockData>("v4", migrateV3ToV4)
  .init(({ params }) => initBlockData(params));

type Ctx = RenderCtx<BlockArgs, BlockData>;

/** The per-replicate frame of the last run. Undefined when the run combined no replicates. */
function replicateCols(ctx: Ctx) {
  return ctx.outputs?.resolve({ field: "replicatePf", allowPermanentAbsence: true })?.getPColumns();
}

function isIntegerValueType(vt: string | undefined): boolean {
  return vt === "Int" || vt === "Long";
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

// Key axes of the datasets this block accepts: bulk and single-cell clonotypes, and variants.
const DATASET_KEY_AXES: ReadonlySet<string> = new Set([
  "pl7.app/vdj/clonotypeKey",
  "pl7.app/vdj/scClonotypeKey",
  "pl7.app/variantKey",
]);

/**
 * The block's workflow input, and the gate on whether it may run at all: a throw here marks the
 * args invalid and disables Run. Named rather than written inline so the staleness contract —
 * which bytes a given state projects — can be asserted in a test.
 */
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

  // Two halves, so `hill_gate_index` can place its settings between them: its args keep
  // v2.4.0's key order as well as its keys.
  const inputs = {
    abundanceRef: data.abundanceRef,
    concentrationColumnRef: data.concentrationColumnRef,
    binColumnRef: data.binColumnRef,
    antigenColumnRef: data.antigenColumnRef,
    sortFractionColumnRef: data.sortFractionColumnRef,
    targetAntigen: data.targetAntigen,
  };
  if (data.replicateColumnRef !== undefined) {
    if (data.replicateFitErrorCutoff !== undefined && data.replicateFitErrorCutoff <= 0)
      throw new Error("Replicate fit-error cutoff must be above 0, or empty to turn it off");
    if (
      data.replicateR2Floor !== undefined &&
      (data.replicateR2Floor < 0 || data.replicateR2Floor > 1)
    )
      throw new Error("Replicate R² floor must be between 0 and 1, or empty to turn it off");
    if (data.replicateCensoredMajority === undefined)
      throw new Error("Replicate censoring majority is required");
    if (data.replicateCensoredMajority < 0.5 || data.replicateCensoredMajority > 1)
      throw new Error("Replicate censoring majority must be between 0.5 and 1");
  }

  const labels: ArgsLabels = {
    // Strings have UI defaults but old stored data may lack them; coerce to
    // empty/safe values rather than throw — they don't break the workflow.
    defaultBlockLabel: data.defaultBlockLabel ?? "Tite-Seq Analysis",
    customBlockLabel: data.customBlockLabel ?? "",
    // The subset's canonical column id, which the workflow stamps as the outputs'
    // `pl7.app/inputSubset`. Absent without a filter, so unfiltered args are unchanged.
    ...(data.filterRef !== undefined && {
      inputFilter: createGlobalPObjectId(data.filterRef.blockId, data.filterRef.name),
    }),
    // Every replicate setting rides on the column: without one the block fits one titration and
    // projects no replicate byte at all.
    ...(data.replicateColumnRef !== undefined && {
      replicateColumnRef: data.replicateColumnRef,
      replicateErrorModel: data.replicateErrorModel,
      ...(data.replicateFitErrorCutoff !== undefined && {
        replicateFitErrorCutoff: data.replicateFitErrorCutoff,
      }),
      ...(data.replicateR2Floor !== undefined && { replicateR2Floor: data.replicateR2Floor }),
      replicateCensoredMajority: data.replicateCensoredMajority,
    }),
  };

  switch (data.equation) {
    case "log_fluorescence":
      return logFluorescenceArgs(data, inputs, labels);
    case "hill_gate_index":
      return hillGateIndexArgs(data, inputs, labels);
    default:
      throw new Error(`Unknown equation: ${data.equation satisfies never}`);
  }
};

function logFluorescenceArgs(data: BlockData, inputs: ArgsInputs, labels: ArgsLabels): BlockArgs {
  // Reject undefined before the range comparison. `undefined <= 0` and `undefined > 1` are both
  // false, so a cleared floor would pass the range check and reach the Python fit.
  if (data.plateauFloor === undefined) throw new Error("Plateau floor is required");
  if (data.plateauFloor <= 0 || data.plateauFloor > 1)
    throw new Error("Plateau floor must be above 0 and at most 1");
  // The log-fluorescence fit's signal is gate fluorescence, so frequency-only data cannot
  // feed it. The Python raises on this; refusing here disables Run and names the reason
  // instead of failing the job after it has been submitted.
  if (data.binColumnRef === undefined)
    throw new Error(
      "The log fluorescence fit requires a FACS bin column. For frequency-only data, choose the Mean bin fit under Fit method.",
    );
  if (data.r2GateEnabled) {
    if (data.r2GateMin === undefined) throw new Error("R² gate threshold is required");
    if (data.r2GateMin < 0 || data.r2GateMin > 1)
      throw new Error("R² gate threshold must be between 0 and 1");
    if (data.kdSeGateMax !== undefined && !(data.kdSeGateMax > 0))
      throw new Error("Maximum Kd,app SE must be above 0, or empty for no limit");
  }
  return {
    ...inputs,
    binFluorescenceColumnRef: data.binFluorescenceColumnRef,
    equation: "log_fluorescence",
    plateauFloor: data.plateauFloor,
    ...(data.r2GateEnabled && { r2GateMin: data.r2GateMin }),
    ...(data.r2GateEnabled && data.kdSeGateMax !== undefined && { kdSeGateMax: data.kdSeGateMax }),
    ...(data.requireAllConcentrations === true && { requireAllConcentrations: true as const }),
    ...labels,
  };
}

function hillGateIndexArgs(data: BlockData, inputs: ArgsInputs, labels: ArgsLabels): BlockArgs {
  // Reject undefined before any range comparison runs. `undefined < 1` returns false, so a
  // stored block that lacks one of these fields would pass validation and send undefined to
  // the Python fit.
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
  if (data.minReadsPerConcentration < 1) throw new Error("Min reads per concentration must be ≥ 1");
  if (data.minConcentrationPoints < 3) throw new Error("Min concentration points must be ≥ 3");
  return {
    ...inputs,
    minReadsPerConcentration: data.minReadsPerConcentration,
    minConcentrationPoints: data.minConcentrationPoints,
    r2ThresholdGood: data.r2ThresholdGood,
    r2ThresholdFailed: data.r2ThresholdFailed,
    nMin: data.nMin,
    nMax: data.nMax,
    // Default the inactive mode's threshold so the args type is satisfied
    // — the workflow only reads the field for the active mode.
    hookEffectThresholdBin:
      data.hookEffectThresholdBin ?? HILL_GATE_INDEX_DEFAULTS.hookEffectThresholdBin,
    hookEffectThresholdNoBin:
      data.hookEffectThresholdNoBin ?? HILL_GATE_INDEX_DEFAULTS.hookEffectThresholdNoBin,
    hookEffectMinReads: data.hookEffectMinReads,
    ...labels,
  };
}

/** The block's pages. Replicates appears only after a run that combined sample replicates. */
export const sectionsFor = (replicates: boolean) => [
  { type: "link" as const, href: "/" as const, label: "Table" },
  ...(replicates
    ? [{ type: "link" as const, href: "/replicates" as const, label: "Replicates" }]
    : []),
  { type: "link" as const, href: "/titration-curves" as const, label: "Titration Curves" },
  { type: "link" as const, href: "/kd-distribution" as const, label: "Kd Distribution" },
  { type: "link" as const, href: "/affinity-vs-fit" as const, label: "Affinity vs Fit Quality" },
];

export const platforma = BlockModelV3.create({ dataModel, kind })

  .templateParams(deriveTemplateParams)

  .args<BlockArgs>(deriveArgs)

  // Datasets, each with its subset columns (e.g. Repertoire Labeling labels) as filters.
  .output(
    "abundanceOptions",
    (ctx): DatasetOption[] =>
      buildDatasetOptions(ctx, {
        primary: (spec) =>
          isPColumnSpec(spec) &&
          spec.annotations?.["pl7.app/isAnchor"] === "true" &&
          spec.axesSpec.length === 2 &&
          spec.axesSpec[0].name === "pl7.app/sampleId" &&
          DATASET_KEY_AXES.has(spec.axesSpec[1].name),
        // Only subsets keyed by the clonotype axis alone: the fit is per clonotype.
        filter: (spec) =>
          isPColumnSpec(spec) &&
          spec.axesSpec.length === 1 &&
          DATASET_KEY_AXES.has(spec.axesSpec[0].name),
        requireEnrichments: false,
      }) ?? [],
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
  // Gate fluorescence is per-sample numeric metadata of the same shape as the sort fraction, so
  // the same filter serves. It is not cross-excluded against the sort fraction: a project may
  // legitimately carry both, and neither constrains the other's values.
  .output("binFluorescenceOptions", (ctx) =>
    filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        isFloatValueType(spec.valueType) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
      { excludeRef: ctx.data.concentrationColumnRef },
    ),
  )
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

  // Which replicate a sample belongs to. A replicate label is a name or a number, so String and
  // integer columns both qualify; a float does not, because a measured quantity is not a label.
  // Cross-excluded against the columns already picked for another role.
  .output("replicateOptions", (ctx) => {
    const picked = [
      ctx.data.concentrationColumnRef,
      ctx.data.binColumnRef,
      ctx.data.sortFractionColumnRef,
      ctx.data.binFluorescenceColumnRef,
      ctx.data.antigenColumnRef,
    ].filter((ref): ref is PlRef => ref !== undefined);
    return filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        (spec.valueType === "String" || isIntegerValueType(spec.valueType)) &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    ).filter((opt) => !picked.some((ref) => plRefsEqual(opt.ref, ref)));
  })

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
    // replicatePf exists only on a run that combined replicates, so its absence is no error.
    const outputKeys = ["summaryPf", "signalPf", "logHandle", "replicatePf"] as const;
    for (const key of outputKeys) {
      const errorAcc = (
        key === "replicatePf"
          ? ctx.outputs?.resolve({ field: key, allowPermanentAbsence: true })
          : ctx.outputs?.resolve(key)
      )?.getError();
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

    // Numeric-field errors surface on each PlNumberField, from its own `required`,
    // `min-value`, `max-value` and `validate` props. The args() throw above gates Run.

    return issues;
  })

  .output("logHandle", (ctx) => ctx.outputs?.resolve("logHandle")?.getLogHandle())

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .outputWithStatus("summaryTable", (ctx) => {
    // One row per variant: the summary columns only. The per-concentration signal, keyed by
    // concentration and replicate as well, would repeat each variant once per measured point; it
    // is plotted on the Titration Curves page and exported in `pf`.
    const summaryCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (summaryCols === undefined) return undefined;

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
      primaryColumns: summaryCols.map((c) => DataColumn.fromColumn(c)),
      columns: resultPoolCols,
      tableState: ctx.data.tableState,
      displayOptions: {
        // Show why each clonotype failed. Results from runs that hid the column still show it.
        visibility: [
          {
            match: { name: [exactly("pl7.app/vdj/fitFailureReason")] },
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

  .output("isEmpty", (ctx) => {
    const pCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols.length === 0;
  })

  .output("binMode", (ctx) => ctx.data.binColumnRef !== undefined)

  // Whether the last run combined sample replicates. The per-replicate page is shown only then.
  .output("replicates", (ctx) => replicateCols(ctx) !== undefined)

  // One row per (clonotype, replicate). Its own page: the extra axis would fan out every row of
  // the main table.
  .outputWithStatus("replicateTable", (ctx) => {
    const cols = replicateCols(ctx);
    if (cols === undefined) return undefined;
    return createPlDataTableV3(ctx, {
      primaryColumns: cols.map((c) => DataColumn.fromColumn(c)),
      columns: undefined,
      tableState: ctx.data.replicateTableState,
      // The per-replicate columns are hidden so downstream pickers offer the combined answer
      // first. This is the page that exists to show them.
      displayOptions: {
        visibility: [
          { match: { name: cols.map((c) => exactly(c.spec.name)) }, visibility: "default" },
        ],
      },
    });
  })

  .outputWithStatus("replicatePfHandle", (ctx): PFrameHandle | undefined => {
    const cols = replicateCols(ctx);
    if (cols === undefined) return undefined;
    return createPFrameForGraphs(ctx, cols);
  })

  .output("replicatePfCols", (ctx) => {
    const cols = replicateCols(ctx);
    if (cols === undefined) return undefined;
    return cols.map((c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec);
  })

  .title(() => "Tite-Seq Analysis")

  .subtitle((ctx) => ctx.data.customBlockLabel || ctx.data.defaultBlockLabel || "")

  .sections((ctx) => sectionsFor(replicateCols(ctx) !== undefined))

  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
