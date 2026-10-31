import type { GraphMakerState } from "@milaboratories/graph-maker";
import type {
  DatasetOption,
  ImportFileHandle,
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlDataTableStateV2,
  PlRef,
} from "@platforma-sdk/model";
import type { BlockParams, Equation } from "@platforma-open/platforma-open.titeseq-analysis.kind";
import { kind } from "@platforma-open/platforma-open.titeseq-analysis.kind";
import type { FitSettings, PresetId } from "./presets";
import { DEFAULT_EQUATION, DEFAULT_PLATEAU_FLOOR, HILL_GATE_INDEX_DEFAULTS } from "./presets";
import type { SettingsImport } from "./settingsFile";
import { readSettingsFile } from "./settingsFile";
import { deriveTemplateParams } from "./templateParams";

export type { BlockParams, Equation } from "@platforma-open/platforma-open.titeseq-analysis.kind";
export { EQUATIONS } from "@platforma-open/platforma-open.titeseq-analysis.kind";
export * from "./presets";
export * from "./settingsFile";
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
  targetAntigen?: string;
};

type BlockLabels = {
  defaultBlockLabel: string;
  customBlockLabel: string;
};

export type BlockArgs = FitInputs &
  BlockLabels &
  Omit<FitSettings, "equation" | "plateauFloor"> & {
    /**
     * Absent means `hill_gate_index` on the wire — deliberately NOT the default a new block
     * takes in `BlockData`, which is `log_fluorescence`. A project that ran before the equation
     * existed projects exactly the args it always projected, so the upgrade does not stale it
     * into re-running a fit that reproduces its own numbers. The workflow and the Python
     * `FitParams` read an absent key the same way.
     */
    equation?: Equation;
    /** Sent only under `log_fluorescence`; no other equation reads a plateau floor. */
    plateauFloor?: number;
    inputFilter?: string;
  };

export type BlockData = FitInputs &
  BlockLabels &
  FitSettings & {
    /** Subset column picked with the abundance in PlDatasetSelector (e.g. a Repertoire Labeling label). */
    filterRef?: PlRef;
    /** The preset last applied, kept so the drawer can name a preset the scientist then edited.
     *  Not authoritative: the settings themselves are. Never reaches the workflow. */
    appliedPreset?: PresetId;
    /** The settings file picked for import. Read through the prerun; applied only on a gesture. */
    importFileHandle?: ImportFileHandle;
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
export type BlockDataV1 = Omit<
  BlockData,
  "equation" | "plateauFloor" | "appliedPreset" | "importFileHandle"
>;

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
  // means the other one — see `BlockArgs.equation`.
  equation: params?.equation ?? DEFAULT_EQUATION,
  plateauFloor: params?.plateauFloor ?? DEFAULT_PLATEAU_FLOOR,
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
 * v1 -> v2. A project created before this release ran the gate-index fit, so it is pinned to it
 * by name. Writing the field is the whole point: an absent `equation` would pick up a new
 * block's default, `log_fluorescence`, and every released project would silently change which
 * equation it fits. The plateau floor comes along at a new block's value and stays unread —
 * `hill_gate_index` has no plateau floor.
 *
 * Nothing else moves. Every stored setting, column reference, label and piece of view state is
 * carried across untouched, so a re-run reproduces the numbers the project already holds.
 */
export const migrateV1ToV2 = (v1: BlockDataV1): BlockData => ({
  ...v1,
  equation: "hill_gate_index",
  plateauFloor: DEFAULT_PLATEAU_FLOOR,
});

export const dataModel = new DataModelBuilder({ kind })
  .from<BlockDataV1>("v1")
  .migrate<BlockData>("v2", migrateV1ToV2)
  .init(({ params }) => initBlockData(params));

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

  // Reject undefined on the numeric tuning params before any range comparisons
  // run. Without these guards, `undefined < 1` (and friends) silently returns
  // false, so a stored block missing one of these fields would pass validation
  // and produce args with undefined values — breaking the Python fit at runtime.
  if (data.equation === undefined) throw new Error("Equation is required");
  if (data.plateauFloor === undefined) throw new Error("Plateau floor is required");
  if (data.plateauFloor <= 0 || data.plateauFloor > 1)
    throw new Error("Plateau floor must be above 0 and at most 1");
  // The log-fluorescence fit's signal is gate fluorescence, so frequency-only data cannot
  // feed it. The Python raises on this; refusing here disables Run and names the reason
  // instead of failing the job after it has been submitted.
  if (data.equation === "log_fluorescence" && data.binColumnRef === undefined)
    throw new Error("The log fluorescence fit requires a FACS bin column");

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
    // Canonicalization, not a shortcut: `hill_gate_index` is what an absent `equation` means
    // on the wire, so a project that stores it projects the bytes v2.4.0 projected and the
    // upgrade does not stale it. The plateau floor rides with the one equation that reads it.
    ...(data.equation !== "hill_gate_index" && {
      equation: data.equation,
      plateauFloor: data.plateauFloor,
    }),
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
    // The subset's canonical column id, which the workflow stamps as the outputs'
    // `pl7.app/inputSubset`. Absent without a filter, so unfiltered args are unchanged.
    ...(data.filterRef !== undefined && {
      inputFilter: createGlobalPObjectId(data.filterRef.blockId, data.filterRef.name),
    }),
  };
};

export const platforma = BlockModelV3.create({ dataModel, kind })

  .templateParams(deriveTemplateParams)

  .args<BlockArgs>(deriveArgs)

  // The settings file is the only thing the prerun carries. Projected on its own so picking one
  // reads it straight away, without staling the block or waiting on an abundance column.
  .prerunArgs((data) => ({ importFileHandle: data.importFileHandle }))

  // isActive forces this output to be evaluated even when nothing subscribes to it. That
  // evaluation is what drives the upload: getImportProgress() causes the underlying context to
  // begin the transfer for each pending ImportFileHandle. Without isActive the output never
  // renders, the upload never starts, and the prerun never resolves — with no error anywhere.
  .output(
    "prerunFileImports",
    (ctx) =>
      Object.fromEntries(
        ctx.prerun
          ?.resolve({ field: "fileImports", assertFieldType: "Input" })
          ?.mapFields((handle, acc) => [handle as ImportFileHandle, acc.getImportProgress()], {
            skipUnresolved: true,
          }) ?? [],
      ),
    { isActive: true },
  )

  // The picked settings file, read and checked here rather than in the UI: the derivation is
  // pure and memoized, and the UI is left with nothing to watch. Applying the settings is a
  // button the scientist presses, never a write this output triggers.
  .output("importedSettings", (ctx) => {
    if (ctx.data.importFileHandle === undefined) return { status: "none" } satisfies SettingsImport;
    const content = ctx.prerun?.traverse({ field: "importedFile" })?.getFileContentAsString();
    if (content === undefined) return { status: "loading" } satisfies SettingsImport;
    return content.map(
      (text): SettingsImport =>
        text === undefined ? { status: "loading" } : readSettingsFile(text),
    );
  })

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

export type BlockOutputs = InferOutputsType<typeof platforma>;
