import type { Equation } from "@platforma-open/platforma-open.titeseq-analysis.kind";

/**
 * Every setting that shapes the fit — the two knobs the equations are chosen and tuned with, and
 * the nine `hill_gate_index` grades its fits by. Inputs, labels and view state are not here.
 *
 * The `hill_gate_index` branch of `BlockArgs` is built on this type, and `BlockData` on
 * `FitSettingsDraft` below, so a setting added here reaches the drawer and the args together. Every value is present.
 */
export type FitSettings = {
  equation: Equation;
  /** The fraction of the gate range a curve must reach to count as saturated. Read by
   *  `log_fluorescence` only. */
  plateauFloor: number;
  /** Whether an in-range fit with R² below `r2GateMin` rates Failed. Read by `log_fluorescence`
   *  only. */
  r2GateEnabled: boolean;
  /** The R² gate's threshold. Read by `log_fluorescence` only, while the gate is on. */
  r2GateMin: number;
  /** The gate's SE(log10 Kd,app) maximum, while the gate is on. Absent sets no limit. Read by
   *  `log_fluorescence` only. */
  kdSeGateMax?: number;
  /** Whether a clonotype with no reads at one of the titration's concentrations, 0 M included,
   *  rates Failed before any fit. Absent is off. Read by both equations. */
  requireAllConcentrations?: boolean;
  minReadsPerConcentration: number;
  minConcentrationPoints: number;
  r2ThresholdGood: number;
  r2ThresholdFailed: number;
  nMin: number;
  nMax: number;
  hookEffectThresholdBin: number;
  hookEffectThresholdNoBin: number;
  hookEffectMinReads: number;
};

/** The settings only `log_fluorescence` reads. */
export type LogFluorescenceSetting = "plateauFloor" | "r2GateEnabled" | "r2GateMin" | "kdSeGateMax";

/** The read-coverage settings both equations read. */
export type ReadCoverageSetting =
  | "minReadsPerConcentration"
  | "minConcentrationPoints"
  | "requireAllConcentrations";

/**
 * The settings as the block stores them, which is what the drawer edits: a number field the
 * scientist cleared reads `undefined`. `BlockData` is built on this type, and `args()` is the
 * single boundary that refuses an incomplete set and so disables Run.
 *
 * `equation` and `r2GateEnabled` stay required: a dropdown and a checkbox have no clear button.
 */
export type FitSettingsDraft = Partial<Omit<FitSettings, "equation" | "r2GateEnabled">> &
  Pick<FitSettings, "equation" | "r2GateEnabled">;

/** The read coverage both equations require, at the values each fit used on its own. */
export const READ_COVERAGE_DEFAULTS = {
  minReadsPerConcentration: 3,
  minConcentrationPoints: 5,
} as const satisfies Pick<FitSettings, "minReadsPerConcentration" | "minConcentrationPoints">;

/** The settings only `hill_gate_index` reads, at the values the block shipped through 2.4.0. */
export const HILL_GATE_INDEX_DEFAULTS = {
  r2ThresholdGood: 0.8,
  r2ThresholdFailed: 0.5,
  nMin: 0.5,
  nMax: 2.0,
  hookEffectThresholdBin: 0.2,
  hookEffectThresholdNoBin: 0.02,
  hookEffectMinReads: 20,
} as const satisfies Omit<FitSettings, "equation" | LogFluorescenceSetting | ReadCoverageSetting>;

export const DEFAULT_EQUATION: Equation = "log_fluorescence";

/** Adams 2016's floor, two thirds of the gate range, to three decimals. */
export const DEFAULT_PLATEAU_FLOOR = 0.667;

/** Off: no rule unverified against a paper's own code sits in the default path. The threshold is
 *  Phillips 2021's. */
export const DEFAULT_R2_GATE = { r2GateEnabled: false, r2GateMin: 0.8 } as const;
