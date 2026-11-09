import type { GraphMakerState } from "@milaboratories/graph-maker";
import type {
  DatasetOption,
  InferOutputsType,
  PColumnIdAndSpec,
  PColumnSpec,
  PFrameHandle,
  PlDataTableFilters,
  PlDataTableStateV2,
  PlRef,
  PTableKey,
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
import type { FitSettings, FitSettingsDraft, LogFluorescenceSetting } from "./fitSettings";
import {
  DEFAULT_EQUATION,
  DEFAULT_PLATEAU_FLOOR,
  DEFAULT_R2_GATE_MIN,
  DEFAULT_KD_SE_GATE_MAX,
  HILL_GATE_INDEX_DEFAULTS,
  READ_COVERAGE_DEFAULTS,
} from "./fitSettings";
import { deriveTemplateParams } from "./templateParams";
import type { InputKind, SynonymousSummary } from "./synonymous";
import {
  halfMatchedParentInput,
  inputKindOf,
  synonymousAaAxis,
  synonymousMessages,
} from "./synonymous";

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
export * from "./fitSettings";
export { deriveTemplateParams };
export type { PTableKey } from "@platforma-sdk/model";
export type { InputKind, SynonymousSummary } from "./synonymous";
export { nounOf } from "./synonymous";
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
  getAxisId,
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
  /**
   * The bin column's values, lowest fluorescence first; position is the gate's rank. A value left
   * out takes no part in the run. Absent reads the column's own integers as the gates, as every
   * block made before the order existed did.
   */
  binOrder?: string[];
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
type ArgsInputs = Omit<
  FitInputs,
  "binFluorescenceColumnRef" | "replicateColumnRef" | "binOrder"
> & {
  /** Each kept bin value's rank, from 1. Present only when the user set an order. */
  binRanks?: Record<string, number>;
};
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
      minReadsPerConcentration: number;
      minConcentrationPoints: number;
      /** Present only when set, so a block with no R² limit names none. */
      r2GateMin?: number;
      /** Present only when set. */
      kdSeGateMax?: number;
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
    /**
     * The variant whose replicates the drill-down shows: the main grid's row key. UI-only, never
     * projected to args. Undefined means the drill-down is closed and its table is not built.
     */
    replicateDrilldownKey?: PTableKey;
    /** Grid state of the drill-down table. UI-only. Optional so no migration is needed. */
    replicateDrilldownTableState?: PlDataTableStateV2;
    graphStateReplicates: GraphMakerState;
    /** Subset column picked with the abundance in PlDatasetSelector (e.g. a Repertoire Labeling label). */
    filterRef?: PlRef;
    tableState: PlDataTableStateV2;
    /** The nucleotide variant table, shown in synonymous mode only. */
    ntTableState: PlDataTableStateV2;
    /** NT Variants page: list only the copies the synonymous graphs list. UI-only, optional. */
    ntSynonymousOnly?: boolean;
    /**
     * The Explorer page, synonymous mode only. UI-only, never projected to args; optional so no
     * migration is needed. The amino-acid row opened, then the nt copy opened inside it.
     */
    explorerTableState?: PlDataTableStateV2;
    explorerAaKey?: PTableKey;
    explorerNtKey?: PTableKey;
    explorerNtTableState?: PlDataTableStateV2;
    explorerReplicateTableState?: PlDataTableStateV2;
    /** The Explorer's combined titration curves, and the affinity classes they show. */
    graphStateExplorerCurves?: GraphMakerState;
    graphStateExplorerScatter?: GraphMakerState;
    explorerCurveClasses?: string[];
    graphStateTitrationCurves: GraphMakerState;
    graphStateKDHistogram: GraphMakerState;
    graphStateAffinityVsFit: GraphMakerState;
    /** Parent mode only. Optional so no migration is needed; the page supplies the default. */
    graphStateVsParent?: GraphMakerState;
    /** Synonymous mode only. Optional so no migration is needed; the page supplies the default. */
    graphStateSynonymousNoise?: GraphMakerState;
    /** The summary pages' graphs. Optional so no migration is needed; each page supplies its default. */
    graphStateReplicateNoise?: GraphMakerState;
    graphStateTitrationEdges?: GraphMakerState;
    /** The Replicates page. Optional so no migration is needed; each view supplies its default. */
    replicateNtTableState?: PlDataTableStateV2;
    graphStateSynonymousCurves?: GraphMakerState;
    graphStateSynonymousKd?: GraphMakerState;
    graphStateReadsByClass?: GraphMakerState;
    settingsOpen: boolean;
  };

/**
 * v1: the block as published through 2.4.0 — one fit, the free-exponent Hill on the read-weighted
 * mean gate index, with no equation to choose. This is the shape every released project stores.
 */
export type BlockDataV1 = Omit<BlockDataV2, "equation" | "plateauFloor">;

/** v2: the equation choice, before the log-fluorescence R² gate. Never released. */
export type BlockDataV2 = Omit<BlockDataV3, "r2GateEnabled" | "r2GateMin">;

/** v6: before the nucleotide variant table, with the error model's 2.5.0 name still allowed. */
export type BlockDataV6 = Omit<BlockData, "ntTableState" | "replicateErrorModel"> & {
  replicateErrorModel: ReplicateErrorModel | LegacyReplicateErrorModel;
};

/** 2.5.0 stored `pop_sd_over_sqrt_n_plus_1` under this name. */
type LegacyReplicateErrorModel = "phillips";

/** v4 and v5: a checkbox turned the R² and SE limits on and off. */
export type BlockDataV5 = BlockDataV6 & { r2GateEnabled: boolean };

/** v3: the R² gate, before sample replicates. Never released. */
export type BlockDataV3 = Omit<
  BlockDataV5,
  | "replicateColumnRef"
  | "replicateErrorModel"
  | "replicateFitErrorCutoff"
  | "replicateR2Floor"
  | "replicateCensoredMajority"
  | "replicateTableState"
  | "graphStateReplicates"
  | "replicateDrilldownKey"
  | "replicateDrilldownTableState"
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

export const PARENT_DELTA_COL = "pl7.app/vdj/kdLog10ParentDelta";

/** The vs Parent histogram of log10 Kd,app vs Parent, on a log count axis. */
export const VS_PARENT_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "bins",
  currentTab: null,
  layersSettings: { bins: { fillColor: "#99E099" } },
  axesSettings: { axisY: { scale: "log" }, other: { binsCount: 30 } },
} as GraphMakerState;

/** The Read Depth box plot: reads in the fit per affinity class, on a log axis. */
export const READ_DEPTH_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "box",
  currentTab: null,
  layersSettings: { box: {} },
  axesSettings: { axisY: { scale: "log" } },
} as GraphMakerState;

/** A histogram on log x and y axes, for the Titration Edges page. */
export const LOG_HISTOGRAM_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "bins",
  currentTab: null,
  layersSettings: { bins: {} },
  axesSettings: { axisX: { scale: "log" }, axisY: { scale: "log" }, other: { binsCount: 30 } },
} as GraphMakerState;

/** The Synonymous Noise and Replicate Noise histograms of a spread in log10 Kd,app. */
/**
 * The synonymous curves: each nt copy's fitted curve at every concentration, a line through the
 * mean over replicates with a standard-error bar.
 */
export const SYNONYMOUS_CURVES_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "line_errorbar",
  currentTab: null,
  layersSettings: { errorbar: { pointsValues: "mean", interval: "se" } },
};

/**
 * The synonymous Kd,app: each nt copy's -log10 Kd,app as the mean over replicates with a
 * standard-error bar.
 */
export const SYNONYMOUS_KD_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "errorbar",
  currentTab: null,
  layersSettings: { errorbar: { pointsValues: "mean", interval: "se" } },
};

/**
 * The Explorer's titration curves: each variant's observed signal and fitted curve at every
 * concentration as a point at the mean over replicates, with a standard-error bar. No lines: a
 * category graph draws the same layers for both series.
 */
export const EXPLORER_CURVES_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "errorbar",
  currentTab: null,
  // Colour follows the series: observed, fitted.
  layersSettings: {
    errorbar: {
      pointsValues: "mean",
      interval: "se",
      lineColor: { type: "secondaryGrouping" },
      dotFill: { type: "secondaryGrouping" },
    },
  },
};

/**
 * The Explorer's titration curves on a log concentration axis: each variant's observed mean over
 * replicates as points and the replicates' mean curve. No error bars: a scatterplot has none.
 */
export const EXPLORER_SCATTER_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "dots",
  currentTab: null,
  layersSettings: { dots: {} },
  axesSettings: { axisX: { scale: "log" } },
};

/** The affinity classes the Explorer's titration curves show until the user picks others. */
export const EXPLORER_CURVE_CLASSES_DEFAULT = ["Good"];

export const NOISE_HISTOGRAM_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "bins",
  currentTab: null,
  layersSettings: { bins: { fillColor: "#99E099" } },
  axesSettings: { other: { binsCount: 30 } },
} as GraphMakerState;

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
  binOrder: params?.binOrder,
  binFluorescenceColumnRef: params?.binFluorescenceColumnRef,
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
  r2GateMin: params?.r2GateMin ?? DEFAULT_R2_GATE_MIN,
  kdSeGateMax: params?.kdSeGateMax ?? DEFAULT_KD_SE_GATE_MAX,
  requireAllConcentrations: params?.requireAllConcentrations,
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
    params?.minReadsPerConcentration ?? READ_COVERAGE_DEFAULTS.minReadsPerConcentration,
  minConcentrationPoints:
    params?.minConcentrationPoints ?? READ_COVERAGE_DEFAULTS.minConcentrationPoints,
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
  ntTableState: createPlDataTableStateV2(),
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
export const migrateV2ToV3 = (v2: BlockDataV2): BlockDataV3 => ({
  ...v2,
  r2GateEnabled: false,
  r2GateMin: 0.8,
});

/**
 * v3 -> v4. Sample replicates arrive unset: no replicate column, both filters off, the default
 * error model and censoring majority, and empty per-replicate view state. The args a v3 project
 * sends do not move, so it does not re-run.
 */
export const migrateV3ToV4 = (v3: BlockDataV3): BlockDataV5 => ({
  ...v3,
  replicateErrorModel: DEFAULT_REPLICATE_ERROR_MODEL,
  replicateCensoredMajority: DEFAULT_CENSORED_MAJORITY,
  replicateTableState: createPlDataTableStateV2(),
  graphStateReplicates: REPLICATES_DEFAULT_STATE,
});

/** The fewest concentration points a fit may require. */
export const MIN_CONCENTRATION_POINTS = 5;

/**
 * v4 -> v5. A stored point count below the minimum is raised to it. A cleared field stays cleared.
 * Raising it changes that project's args, so the block goes stale.
 */
export const migrateV4ToV5 = (v4: BlockDataV5): BlockDataV5 => ({
  ...v4,
  minConcentrationPoints:
    v4.minConcentrationPoints === undefined
      ? undefined
      : Math.max(v4.minConcentrationPoints, MIN_CONCENTRATION_POINTS),
});

/**
 * v5 -> v6. The checkbox is gone: a limit applies when it is set. A project with the checkbox off
 * clears both limits, and one with it on keeps them, so the args do not move and it does not
 * re-run.
 */
export const migrateV5ToV6 = ({ r2GateEnabled, ...v5 }: BlockDataV5): BlockDataV6 => ({
  ...v5,
  r2GateMin: r2GateEnabled ? v5.r2GateMin : undefined,
  kdSeGateMax: r2GateEnabled ? v5.kdSeGateMax : undefined,
});

/** v6 -> v7. Adds the nucleotide variant table's state, empty, and renames the stored error model. */
export const migrateV6ToV7 = (v6: BlockDataV6): BlockData => ({
  ...v6,
  replicateErrorModel:
    v6.replicateErrorModel === "phillips" ? "pop_sd_over_sqrt_n_plus_1" : v6.replicateErrorModel,
  ntTableState: createPlDataTableStateV2(),
});

export const dataModel = new DataModelBuilder({ kind })
  .from<BlockDataV1>("v1")
  .migrate<BlockDataV2>("v2", migrateV1ToV2)
  .migrate<BlockDataV3>("v3", migrateV2ToV3)
  .migrate<BlockDataV5>("v4", migrateV3ToV4)
  .migrate<BlockDataV5>("v5", migrateV4ToV5)
  .migrate<BlockDataV6>("v6", migrateV5ToV6)
  .migrate<BlockData>("v7", migrateV6ToV7)
  .init(({ params }) => initBlockData(params));

type Ctx = RenderCtx<BlockArgs, BlockData>;

/** The per-replicate frame of the last run. Undefined when the run combined no replicates. */
function replicateCols(ctx: Ctx) {
  return ctx.outputs?.resolve({ field: "replicatePf", allowPermanentAbsence: true })?.getPColumns();
}

const REPLICATE_AXIS = "pl7.app/vdj/replicate";
/** Each nt copy's amino-acid variant, on the nt axis: the synonymous graphs tab by it. */
const AA_VARIANT_ID = "pl7.app/vdj/aaVariantId";
/** Each nt copy's amino-acid variant key, on the nt axis: the Explorer filters a variant's copies by it. */
const AA_VARIANT_KEY = "pl7.app/vdj/aaVariantKey";

/**
 * The filter that narrows the per-replicate frame to one variant. The variant axis is whichever
 * axis is not the replicate one, taken from the frame's own spec so the domain matches exactly.
 * Undefined when nothing is selected or the frame has no such axis.
 */
export function replicateDrilldownFilter(
  specs: PColumnSpec[],
  key: PTableKey | undefined,
): PlDataTableFilters | undefined {
  const value = key?.[0];
  if (value === undefined || value === null) return undefined;
  const axis = specs
    .find((s) => s.name === "pl7.app/vdj/kd")
    ?.axesSpec.find((a) => a.name !== REPLICATE_AXIS);
  if (axis === undefined) return undefined;
  return {
    type: "and",
    filters: [
      {
        type: "patternEquals",
        column: { type: "axis", id: getAxisId(axis) },
        value: String(value),
      },
    ],
  };
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

/** A stable string key for a ref, for maps keyed by column. */
export const refKey = (ref: PlRef): string => `${ref.blockId}/${ref.name}`;

/** Each bin value's rank: its position in the order, from 1. */
export const binRanks = (order: string[]): Record<string, number> =>
  Object.fromEntries(order.map((value, i) => [value, i + 1]));

/** A column's values in the order a new pick seeds: numerically when all are numbers. */
export const defaultBinOrder = (values: string[]): string[] => {
  const unique = [...new Set(values)];
  const numeric = unique.every((v) => v.trim() !== "" && Number.isFinite(Number(v)));
  return numeric
    ? unique.sort((a, b) => Number(a) - Number(b))
    : unique.sort((a, b) => a.localeCompare(b));
};

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
  if (data.binColumnRef !== undefined && data.binOrder !== undefined && data.binOrder.length === 0)
    throw new Error("Keep at least one FACS bin in the bin order");

  // Two halves, so `hill_gate_index` can place its settings between them: its args keep
  // v2.4.0's key order as well as its keys.
  const inputs = {
    abundanceRef: data.abundanceRef,
    concentrationColumnRef: data.concentrationColumnRef,
    binColumnRef: data.binColumnRef,
    antigenColumnRef: data.antigenColumnRef,
    sortFractionColumnRef: data.sortFractionColumnRef,
    targetAntigen: data.targetAntigen,
    // Absent without an order, so a block that never set one projects the args it always did.
    ...(data.binColumnRef !== undefined &&
      data.binOrder !== undefined && { binRanks: binRanks(data.binOrder) }),
  };
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
    // projects no replicate byte at all. While the Sample replicates section is hidden, every
    // setting in it is pinned to its default — sample SD / sqrt(n), both filters off, a simple
    // majority — so a value a block stored before the section went away is not used unseen.
    // Restore the `data.replicate*` fields, and their checks, when that section comes back.
    ...(data.replicateColumnRef !== undefined && {
      replicateColumnRef: data.replicateColumnRef,
      replicateErrorModel: DEFAULT_REPLICATE_ERROR_MODEL,
      replicateCensoredMajority: DEFAULT_CENSORED_MAJORITY,
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

/**
 * Read coverage, which both equations read: a point below the read floor is dropped, and a
 * clonotype left with fewer than the point count is refused before any fit. Narrowed here because
 * a check in the caller does not narrow `data` inside each equation's builder.
 */
function readCoverage(data: BlockData): {
  minReadsPerConcentration: number;
  minConcentrationPoints: number;
} {
  if (data.minReadsPerConcentration === undefined)
    throw new Error("Min reads per concentration is required");
  if (data.minConcentrationPoints === undefined)
    throw new Error("Min concentration points is required");
  if (data.minReadsPerConcentration < 1) throw new Error("Min reads per concentration must be ≥ 1");
  if (data.minConcentrationPoints < MIN_CONCENTRATION_POINTS)
    throw new Error(`Min concentration points must be ≥ ${MIN_CONCENTRATION_POINTS}`);
  return {
    minReadsPerConcentration: data.minReadsPerConcentration,
    minConcentrationPoints: data.minConcentrationPoints,
    // `requireAllConcentrations` is withheld: absent fits a titration with a gap.
  };
}

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
      "The log fluorescence fit requires a FACS bin column. For frequency-only data, choose the Mean bin fit in Fit method equation.",
    );
  if (data.r2GateMin !== undefined && !(data.r2GateMin >= 0 && data.r2GateMin <= 1))
    throw new Error("Minimum R² must be between 0 and 1, or empty for no limit");
  if (data.kdSeGateMax !== undefined && !(data.kdSeGateMax > 0))
    throw new Error("Maximum Kd,app SE must be above 0, or empty for no limit");
  return {
    ...inputs,
    binFluorescenceColumnRef: data.binFluorescenceColumnRef,
    equation: "log_fluorescence",
    plateauFloor: data.plateauFloor,
    ...readCoverage(data),
    ...(data.r2GateMin !== undefined && { r2GateMin: data.r2GateMin }),
    ...(data.kdSeGateMax !== undefined && { kdSeGateMax: data.kdSeGateMax }),
    ...labels,
  };
}

function hillGateIndexArgs(data: BlockData, inputs: ArgsInputs, labels: ArgsLabels): BlockArgs {
  // Reject undefined before any range comparison runs. `undefined < 1` returns false, so a
  // stored block that lacks one of these fields would pass validation and send undefined to
  // the Python fit. Read coverage is checked once, before the equation is chosen.
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
  return {
    ...inputs,
    ...readCoverage(data),
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

/** The block's pages. Per-replicate values are read by clicking a variant in the main table. */
/** The amino-acid frame of the last run; undefined outside synonymous mode. */
function variantSummaryCols(ctx: Ctx) {
  return ctx.outputs
    ?.resolve({ field: "variantSummaryPf", allowPermanentAbsence: true })
    ?.getPColumns();
}

function synonymousSummary(ctx: Ctx): SynonymousSummary | undefined {
  // A summary from a run on another dataset says nothing about this one.
  if (!plRefsEqualOrBothAbsent(ctx.activeArgs?.abundanceRef, ctx.data.abundanceRef))
    return undefined;
  return ctx.outputs
    ?.resolve({ field: "synonymousSummary", allowPermanentAbsence: true })
    ?.getDataAsJson<SynonymousSummary>();
}

/** What the last run must tell the user. Empty for a run on another dataset, or before any run. */
function runWarnings(ctx: Ctx): string[] {
  if (!plRefsEqualOrBothAbsent(ctx.activeArgs?.abundanceRef, ctx.data.abundanceRef)) return [];
  return (
    ctx.outputs
      ?.resolve({ field: "runWarnings", allowPermanentAbsence: true })
      ?.getDataAsJson<string[]>() ?? []
  );
}

function plRefsEqualOrBothAbsent(a: PlRef | undefined, b: PlRef | undefined): boolean {
  if (a === undefined || b === undefined) return a === b;
  return plRefsEqual(a, b);
}

/** The picked abundance's amino-acid axis when synonymous mode will turn on, read from the pool. */
function poolSynonymousAaAxis(ctx: Ctx) {
  if (ctx.data.abundanceRef === undefined) return undefined;
  const ntAxis = ctx.resultPool.getPColumnSpecByRef(ctx.data.abundanceRef)?.axesSpec[1];
  if (ntAxis === undefined) return undefined;
  const specs = ctx.resultPool.getSpecs().entries.map((e) => e.obj);
  const aaAxis = synonymousAaAxis(specs, ntAxis);
  return aaAxis === undefined ? undefined : { specs, aaAxis };
}

/** Read from the result pool, so the warning shows before the first run. */
function halfMatchedParent(ctx: Ctx) {
  const found = poolSynonymousAaAxis(ctx);
  return found === undefined ? undefined : halfMatchedParentInput(found.specs, found.aaAxis);
}

/**
 * The block's pages. The nucleotide variant table exists in synonymous mode only, beside the
 * amino-acid table. Each entry stays an `as const` literal: the routes are typed from the hrefs.
 */
/**
 * The log-fluorescence fit places each gate at its measured log10 fluorescence. With no gate
 * fluorescence column it places the gates 0.65 decades apart, and the run gives no other sign of it.
 */
export const GATE_FLUORESCENCE_MISSING =
  "Gate fluorescence is not set. The fit assumes that the gates are 0.65 decades apart. " +
  "If they are not, Kd,app can be offset. To use the measured gate positions, select a " +
  "Gate fluorescence column in the settings.";

/** Advice that the settings alone cause, before any run. The settings drawer shows it. */
export const settingsWarnings = (data: BlockData): string[] =>
  data.equation === "log_fluorescence" &&
  data.binColumnRef !== undefined &&
  data.binFluorescenceColumnRef === undefined
    ? [GATE_FLUORESCENCE_MISSING]
    : [];

/**
 * What a finished run found. Undefined before the first run. `synonymousCopies`: some amino-acid
 * variant has two or more nucleotide copies.
 */
export type RunShape = {
  synonymous: boolean;
  synonymousCopies?: boolean;
  parentMode?: boolean;
  replicates?: boolean;
};

/**
 * The QC page's views, one tab each, in this order. A view that reads data the run did not produce
 * gets no tab: the replicate spread needs replicates, the parent delta parent mode, the synonymous
 * spread an amino-acid variant with two or more nucleotide copies.
 */
export const qcViewsFor = (run: RunShape) => {
  const { synonymous, synonymousCopies = false, parentMode = false, replicates = false } = run;
  return [
    { value: "read-depth" as const, label: "Read Depth" },
    { value: "affinity-vs-fit" as const, label: "Affinity vs Fit Quality" },
    { value: "titration-curves" as const, label: "Titration Curves" },
    { value: "titration-edges" as const, label: "Titration Edges" },
    ...(replicates ? [{ value: "replicate-noise" as const, label: "Replicate Noise" }] : []),
    ...(parentMode ? [{ value: "vs-parent" as const, label: "Kd,app vs Parent" }] : []),
    ...(synonymous && synonymousCopies
      ? [{ value: "synonymous-noise" as const, label: "Synonymous Noise" }]
      : []),
    { value: "kd-distribution" as const, label: "Kd Distribution" },
  ];
};

export type QcView = ReturnType<typeof qcViewsFor>[number]["value"];

/**
 * The Replicates page's views, one tab each: every variant's value in every replicate, then in
 * synonymous mode each nt copy's, then the copies' curves where some variant has two.
 */
export const replicateViewsFor = (run: RunShape) => {
  const { synonymous, synonymousCopies = false } = run;
  return [
    {
      value: "variant-replicates" as const,
      label: synonymous ? "AA by Replicate" : "By Replicate",
    },
    ...(synonymous ? [{ value: "nt-replicates" as const, label: "NT by Replicate" }] : []),
    ...(synonymous && synonymousCopies
      ? [
          { value: "synonymous-curves" as const, label: "Synonymous Curves" },
          { value: "synonymous-kd" as const, label: "Synonymous Kd,app" },
        ]
      : []),
  ];
};

export type ReplicateView = ReturnType<typeof replicateViewsFor>[number]["value"];

/**
 * The pages: the two tables, the Replicates page where the run had replicates, then one QC page.
 * Before the first run only Main shows: every other page reads the run's results.
 */
export const sectionsFor = (run: RunShape | undefined) => {
  const main = { type: "link" as const, href: "/" as const, label: "Main" };
  if (run === undefined) return [main];
  return [
    main,
    ...(run.synonymous
      ? [
          { type: "link" as const, href: "/nt-variants" as const, label: "NT Variants" },
          { type: "link" as const, href: "/explorer" as const, label: "Explorer" },
        ]
      : []),
    ...(run.replicates
      ? [{ type: "link" as const, href: "/replicates" as const, label: "Replicates" }]
      : []),
    { type: "link" as const, href: "/qc" as const, label: "QC" },
  ];
};

/** The shape of the last run, read from its outputs; undefined before the first run. */
function runShape(ctx: Ctx): RunShape | undefined {
  if (ctx.outputs?.resolve({ field: "summaryPf", allowPermanentAbsence: true }) === undefined)
    return undefined;
  const cols = variantSummaryCols(ctx);
  // Read unfiltered: the run's shape holds whatever dataset it ran on. Every mapped copy belongs to
  // one amino-acid variant, so more mapped copies than variants means some variant has two.
  const summary = ctx.outputs
    ?.resolve({ field: "synonymousSummary", allowPermanentAbsence: true })
    ?.getDataAsJson<SynonymousSummary>();
  return {
    synonymous: cols !== undefined,
    synonymousCopies:
      summary !== undefined &&
      summary.ntVariants - summary.unmappedNtVariants - summary.ambiguousNtVariants >
        summary.aaVariants,
    parentMode: cols?.some((c) => c.spec.name === PARENT_DELTA_COL) ?? false,
    replicates: replicateCols(ctx) !== undefined,
  };
}

type FitColumns = NonNullable<ReturnType<typeof variantSummaryCols>>;

/**
 * A fit results table: `summaryCols` as primary columns, joined with the result pool's columns on
 * the Kd,app column's axis.
 */
/**
 * The amino-acid variant's ID beside each nucleotide variant: the amino-acid axis label, reached
 * through the nt -> aa linker. Shown second, after the nucleotide variant's own ID.
 */
function aaVariantIdColumns(ctx: Ctx, kdSpec: PColumnSpec) {
  return ColumnsCollection(["result_pool"], { ctx: ctx.ctx })
    .discover({
      anchors: { main: kdSpec },
      mode: "enrichment",
      maxHops: 1,
      include: {
        name: [exactly("pl7.app/label")],
        axes: [
          {
            name: [exactly("pl7.app/variantKey")],
            domain: { "pl7.app/alphabet": [exactly("aminoacid")] },
          },
        ],
      },
    })
    .getColumns()
    .slice(0, 1)
    .map((c) =>
      c.withSpecs({
        annotations: {
          "pl7.app/label": "AA Variant Id",
          "pl7.app/table/visibility": "default",
          "pl7.app/table/orderPriority": "99000",
        },
      }),
    );
}

function fitTable(
  ctx: Ctx,
  summaryCols: FitColumns,
  tableState: PlDataTableStateV2 | undefined,
  {
    linkAaVariant = false,
    filters,
  }: { linkAaVariant?: boolean; filters?: PlDataTableFilters } = {},
) {
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
    .map((c) => {
      // Behind the block's own columns (176 and up), in the producer's order.
      const priority = Number(c.getSpec().annotations?.["pl7.app/table/orderPriority"]);
      return c.withSpecs({
        annotations: {
          "pl7.app/table/visibility": "optional",
          ...(Number.isFinite(priority) && {
            "pl7.app/table/orderPriority": String(priority / 10_000),
          }),
        },
      });
    });

  return createPlDataTableV3(ctx, {
    primaryColumns: summaryCols.map((c) => DataColumn.fromColumn(c)),
    columns: linkAaVariant
      ? [...aaVariantIdColumns(ctx, kdCol.spec), ...resultPoolCols]
      : resultPoolCols,
    tableState,
    filters,
  });
}

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
        (isIntegerValueType(spec.valueType) || spec.valueType === "String") &&
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

  // Each bin option's distinct values, keyed by the option's ref, so a pick can seed the order in
  // the same gesture: an output on the picked column alone would arrive a round trip late.
  .output("binValuesByRef", (ctx): Record<string, string[]> | undefined => {
    const options = filterPopulatedOptions(
      ctx,
      (spec) =>
        isPColumnSpec(spec) &&
        (isIntegerValueType(spec.valueType) || spec.valueType === "String") &&
        spec.axesSpec.length === 1 &&
        spec.axesSpec[0].name === "pl7.app/sampleId",
    );
    const out: Record<string, string[]> = {};
    for (const option of options) {
      const values = ctx.resultPool
        .getDataByRef(option.ref)
        ?.data?.getDataAsJson<Record<string, Record<string, string | number | null>>>()?.["data"];
      if (values === undefined) continue;
      out[refKey(option.ref)] = defaultBinOrder(
        Object.values(values)
          .filter((v): v is string | number => v !== null && v !== "")
          .map(String),
      );
    }
    return out;
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
    // replicatePf exists only on a run that combined replicates, the next three in synonymous mode
    // only, replicateNtPf only in both, and runWarnings and concentrationOrder not on a run made
    // before they existed, so their absence is no error.
    const optionalKeys = new Set([
      "replicatePf",
      "variantSummaryPf",
      "synonymousSummary",
      "ntVariantPf",
      "replicateNtPf",
      "runWarnings",
      "concentrationOrder",
    ]);
    const outputKeys = [
      "summaryPf",
      "signalPf",
      "logHandle",
      "replicatePf",
      "variantSummaryPf",
      "synonymousSummary",
      "ntVariantPf",
      "replicateNtPf",
      "runWarnings",
      "concentrationOrder",
    ] as const;
    for (const key of outputKeys) {
      const errorAcc = (
        optionalKeys.has(key)
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

    for (const message of synonymousMessages(halfMatchedParent(ctx), synonymousSummary(ctx)))
      issues.push({ severity: "warning", message });

    for (const message of runWarnings(ctx)) issues.push({ severity: "warning", message });

    return issues;
  })

  // Advice about the settings alone, shown in the settings drawer only, not on the results pages.
  .output("settingsWarnings", (ctx): string[] => settingsWarnings(ctx.data))

  .output("logHandle", (ctx) => ctx.outputs?.resolve("logHandle")?.getLogHandle())

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .output("inputKind", (ctx): InputKind | undefined => {
    if (ctx.data.abundanceRef === undefined) return undefined;
    const axis = ctx.resultPool.getPColumnSpecByRef(ctx.data.abundanceRef)?.axesSpec[1];
    return axis === undefined ? undefined : inputKindOf(axis);
  })

  // In synonymous mode the table shows one row per amino-acid variant. Each nucleotide variant's own
  // fit stays in summaryPf, which the graphs read.
  .output("synonymous", (ctx) => variantSummaryCols(ctx) !== undefined)

  .outputWithStatus("summaryTable", (ctx) => {
    // One row per variant: the summary columns only. The per-concentration signal, keyed by
    // concentration and replicate as well, would repeat each variant once per measured point; it
    // is plotted on the Titration Curves page and exported in `pf`. In synonymous mode the rows are
    // amino-acid variants.
    const summaryCols = variantSummaryCols(ctx) ?? ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (summaryCols === undefined) return undefined;

    return fitTable(ctx, summaryCols, ctx.data.tableState);
  })

  // Synonymous mode only: one row per nucleotide variant, its own fit and its nt-grain columns.
  .outputWithStatus("ntVariantTable", (ctx) => {
    if (variantSummaryCols(ctx) === undefined) return undefined;
    const summaryCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (summaryCols === undefined) return undefined;
    const ntCols =
      ctx.outputs?.resolve({ field: "ntVariantPf", allowPermanentAbsence: true })?.getPColumns() ??
      [];
    const shown = ntCols.find((c) => c.spec.name === AA_VARIANT_ID);
    return fitTable(ctx, [...summaryCols, ...ntCols], ctx.data.ntTableState, {
      linkAaVariant: true,
      filters:
        ctx.data.ntSynonymousOnly && shown
          ? {
              type: "and",
              filters: [{ type: "isNotNA", column: { type: "column", id: shown.id } }],
            }
          : undefined,
    });
  })

  // The Explorer page, synonymous mode only: the amino-acid variants, then one variant's nt
  // copies, then one copy's replicates.
  .outputWithStatus("explorerTable", (ctx) => {
    const cols = variantSummaryCols(ctx);
    if (cols === undefined) return undefined;
    return fitTable(ctx, cols, ctx.data.explorerTableState);
  })

  // No table while no amino-acid variant is open: unfiltered it would be every copy.
  .outputWithStatus("explorerNtTable", (ctx) => {
    const aa = ctx.data.explorerAaKey?.[0];
    if (aa === undefined || aa === null) return undefined;
    const summaryCols = ctx.outputs?.resolve("summaryPf")?.getPColumns();
    const ntCols = ctx.outputs
      ?.resolve({ field: "ntVariantPf", allowPermanentAbsence: true })
      ?.getPColumns();
    const keyCol = ntCols?.find((c) => c.spec.name === AA_VARIANT_KEY);
    if (summaryCols === undefined || ntCols === undefined || keyCol === undefined) return undefined;
    // No linker join: every row is the same amino-acid variant, and the row key stays the nt key.
    return fitTable(ctx, [...summaryCols, ...ntCols], ctx.data.explorerNtTableState, {
      filters: {
        type: "and",
        filters: [
          {
            type: "patternEquals",
            column: { type: "column", id: keyCol.id },
            value: String(aa),
          },
        ],
      },
    });
  })

  // One nt copy's replicates. No table while no copy is open, or without replicates.
  .outputWithStatus("explorerReplicateTable", (ctx) => {
    const cols = ctx.outputs
      ?.resolve({ field: "replicateNtPf", allowPermanentAbsence: true })
      ?.getPColumns()
      ?.filter((c) => c.spec.name !== AA_VARIANT_ID && c.spec.name !== AA_VARIANT_KEY);
    if (cols === undefined) return undefined;
    const filters = replicateDrilldownFilter(
      cols.map((c) => c.spec),
      ctx.data.explorerNtKey,
    );
    if (filters === undefined) return undefined;
    return createPlDataTableV3(ctx, {
      primaryColumns: cols.map((c) => DataColumn.fromColumn(c)),
      columns: undefined,
      tableState: ctx.data.explorerReplicateTableState,
      displayOptions: {
        visibility: [
          { match: { name: [exactly("pl7.app/label")] }, visibility: "hidden" },
          { match: { name: cols.map((c) => exactly(c.spec.name)) }, visibility: "default" },
        ],
      },
      filters,
    });
  })

  .outputWithStatus("titrationCurvesPf", (ctx): PFrameHandle | undefined => {
    const signalCols = ctx.outputs?.resolve("signalPf")?.getPColumns();
    if (signalCols === undefined) return undefined;
    const summaryCols = ctx.outputs?.resolve("summaryPf")?.getPColumns() ?? [];
    return createPFrameForGraphs(ctx, [...signalCols, ...summaryCols]);
  })

  // The non-zero concentrations as the String axis holds them, lowest first: the order the
  // category graphs show them in. Undefined on a run made before it existed.
  .output("concentrationOrder", (ctx) =>
    ctx.outputs
      ?.resolve({ field: "concentrationOrder", allowPermanentAbsence: true })
      ?.getDataAsJson<string[]>(),
  )

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

  // What the last run produced, which decides the QC page's tabs. Undefined before the first run.
  .output("runShape", (ctx) => runShape(ctx))

  // One row per (clonotype, replicate). Its own page: the extra axis would fan out every row of
  // the main table.
  .outputWithStatus("replicateTable", (ctx) => {
    const cols = replicateCols(ctx);
    if (cols === undefined) return undefined;
    // The frame's AA Variant Id is for the graphs; the table reaches the label through the linker.
    const tableCols = cols.filter(
      (c) => c.spec.name !== AA_VARIANT_ID && c.spec.name !== AA_VARIANT_KEY,
    );
    return createPlDataTableV3(ctx, {
      primaryColumns: tableCols.map((c) => DataColumn.fromColumn(c)),
      columns: undefined,
      tableState: ctx.data.replicateTableState,
      // The per-replicate columns are hidden so downstream pickers offer the combined answer
      // first. This is the page that exists to show them.
      displayOptions: {
        visibility: [
          { match: { name: tableCols.map((c) => exactly(c.spec.name)) }, visibility: "default" },
        ],
      },
    });
  })

  // The main grid's key axis. Feeds `showCellButtonForAxisId`, which renders the button on every
  // displayed column that sits on this axis alone. Both columns carry the same axis, so the id is
  // the same either way.
  .output("summaryKeyAxisId", (ctx) => {
    // In synonymous mode the grid's rows are amino-acid variants, and so is the per-replicate frame.
    const cols = variantSummaryCols(ctx) ?? ctx.outputs?.resolve("summaryPf")?.getPColumns();
    const axis = (
      cols?.find((c) => c.spec.name === "pl7.app/vdj/affinityClass") ??
      cols?.find((c) => c.spec.name === "pl7.app/vdj/kd")
    )?.spec.axesSpec[0];
    return axis === undefined ? undefined : getAxisId(axis);
  })

  // No `replicateDrilldownLabel`: the model cannot name the variant. The clonotype axis value is
  // an opaque key, its label column is Parquet-partitioned so `resultPool.findLabels` throws on
  // it, and neither table click event carries the row's cells. The UI reads the label from this
  // table's hidden label column instead.

  // One variant's replicates, opened from the main table. No table while nothing is selected:
  // unfiltered it would be every variant's replicates.
  .outputWithStatus("replicateDrilldownTable", (ctx) => {
    const cols = replicateCols(ctx);
    if (cols === undefined) return undefined;
    const filters = replicateDrilldownFilter(
      cols.map((c) => c.spec),
      ctx.data.replicateDrilldownKey,
    );
    if (filters === undefined) return undefined;
    // Every per-replicate column is shown. The variant's label column is hidden: the panel title
    // names the variant, and the table holds one variant only.
    return createPlDataTableV3(ctx, {
      primaryColumns: cols.map((c) => DataColumn.fromColumn(c)),
      columns: undefined,
      tableState: ctx.data.replicateDrilldownTableState,
      displayOptions: {
        visibility: [
          { match: { name: [exactly("pl7.app/label")] }, visibility: "hidden" },
          { match: { name: cols.map((c) => exactly(c.spec.name)) }, visibility: "default" },
        ],
      },
      filters,
    });
  })

  .outputWithStatus("replicatePfHandle", (ctx): PFrameHandle | undefined => {
    const cols = replicateCols(ctx);
    if (cols === undefined) return undefined;
    return createPFrameForGraphs(ctx, cols);
  })

  // The reported grain alone: amino-acid variants in synonymous mode, the input's variants
  // otherwise. The summary pages' graphs read it.
  .outputWithStatus("reportedPfHandle", (ctx): PFrameHandle | undefined => {
    const cols = variantSummaryCols(ctx) ?? ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (cols === undefined) return undefined;
    return createPFrameForGraphs(ctx, cols);
  })

  .output("reportedPfCols", (ctx) => {
    const cols = variantSummaryCols(ctx) ?? ctx.outputs?.resolve("summaryPf")?.getPColumns();
    if (cols === undefined) return undefined;
    return cols.map((c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec);
  })

  .output("synonymousSummary", (ctx) => synonymousSummary(ctx))

  // Synonymous mode with replicates: one row per (nt copy, replicate), its own fit, beside the
  // amino-acid variant it belongs to.
  .outputWithStatus("replicateNtTable", (ctx) => {
    const cols = ctx.outputs
      ?.resolve({ field: "replicateNtPf", allowPermanentAbsence: true })
      ?.getPColumns();
    if (cols === undefined) return undefined;
    const kdCol = cols.find((c) => c.spec.name === "pl7.app/vdj/kd");
    const ntAxis = kdCol?.spec.axesSpec.find((a) => a.name !== REPLICATE_AXIS);
    if (kdCol === undefined || ntAxis === undefined) return undefined;
    return createPlDataTableV3(ctx, {
      primaryColumns: cols.map((c) => DataColumn.fromColumn(c)),
      // The linker reaches the amino-acid label from the nt axis alone: an anchor that also
      // carries the replicate axis finds nothing.
      columns: aaVariantIdColumns(ctx, { ...kdCol.spec, axesSpec: [ntAxis] }),
      tableState: ctx.data.replicateNtTableState,
      displayOptions: {
        visibility: [
          { match: { name: cols.map((c) => exactly(c.spec.name)) }, visibility: "default" },
        ],
      },
    });
  })

  // Synonymous mode with replicates: the per-replicate nt fits and each copy's amino-acid
  // variant, for the Synonymous Kd,app graph.
  .outputWithStatus("synonymousKdPf", (ctx): PFrameHandle | undefined => {
    const cols = ctx.outputs
      ?.resolve({ field: "replicateNtPf", allowPermanentAbsence: true })
      ?.getPColumns();
    if (cols === undefined) return undefined;
    return createPFrameForGraphs(ctx, cols);
  })

  .output("synonymousKdPfCols", (ctx) =>
    ctx.outputs
      ?.resolve({ field: "replicateNtPf", allowPermanentAbsence: true })
      ?.getPColumns()
      ?.map((c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec),
  )

  .output("replicatePfCols", (ctx) => {
    const cols = replicateCols(ctx);
    if (cols === undefined) return undefined;
    return cols.map((c) => ({ columnId: c.id, spec: c.spec }) satisfies PColumnIdAndSpec);
  })

  .title(() => "Tite-Seq Analysis")

  .subtitle((ctx) => ctx.data.customBlockLabel || ctx.data.defaultBlockLabel || "")

  .sections((ctx) => sectionsFor(runShape(ctx)))

  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
