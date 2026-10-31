import type { Equation, FitSettings } from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  EQUATIONS,
  HILL_GATE_INDEX_DEFAULTS,
} from "@platforma-open/platforma-open.titeseq-analysis.model";

// Exhaustive over the model's equation union: a new equation fails to compile until it is named.
const EQUATION_LABELS: Record<Equation, string> = {
  log_fluorescence: "Log fluorescence fit",
  hill_gate_index: "Mean gate index fit (previous default)",
};

export const EQUATION_OPTIONS: readonly { readonly value: Equation; readonly label: string }[] =
  EQUATIONS.map((value) => ({ value, label: EQUATION_LABELS[value] }));

export function equationLabel(equation: Equation): string {
  return EQUATION_LABELS[equation];
}

/** Field names as the settings drawer spells them, for the import preview and the fit summary. */
export const SETTINGS_FIELD_LABELS: Record<keyof FitSettings, string> = {
  equation: "Equation",
  plateauFloor: "Plateau floor",
  minReadsPerConcentration: "Min reads per concentration",
  minConcentrationPoints: "Min concentration points",
  r2ThresholdGood: "R² threshold — Good",
  r2ThresholdFailed: "R² threshold — Failed",
  nMin: "Hill coefficient — min",
  nMax: "Hill coefficient — max",
  hookEffectThresholdBin: "Signal drop threshold (bin mode)",
  hookEffectThresholdNoBin: "Signal drop threshold (frequency mode)",
  hookEffectMinReads: "Min reads for hook check",
};

/** Every setting but the equation; all of them numeric. */
export type NumberSetting = Exclude<keyof FitSettings, "equation">;

// Keys come from SETTINGS_FIELD_LABELS itself, which is exhaustive over FitSettings.
export const NUMBER_SETTINGS: readonly NumberSetting[] = Object.keys(SETTINGS_FIELD_LABELS).filter(
  (key): key is NumberSetting => key !== "equation" && key in SETTINGS_FIELD_LABELS,
);

type HillOnlyField = keyof typeof HILL_GATE_INDEX_DEFAULTS;

// Keys come from HILL_GATE_INDEX_DEFAULTS itself, which is the model's own statement of which
// settings only `hill_gate_index` reads.
const HILL_ONLY_FIELDS: readonly HillOnlyField[] = Object.keys(HILL_GATE_INDEX_DEFAULTS).filter(
  (key): key is HillOnlyField => key in HILL_GATE_INDEX_DEFAULTS,
);

/**
 * The settings each equation reads. The rest stay in `BlockData` untouched: the drawer hides them
 * and the fit summary leaves them out, because the running equation never consults them.
 */
export const FIELDS_READ_BY: Record<Equation, readonly (keyof FitSettings)[]> = {
  log_fluorescence: ["plateauFloor"],
  hill_gate_index: HILL_ONLY_FIELDS,
};
