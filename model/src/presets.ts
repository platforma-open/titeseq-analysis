import type { Equation } from "@platforma-open/platforma-open.titeseq-analysis.kind";

/**
 * Every setting that shapes the fit — the two knobs the equations are chosen and tuned with, and
 * the nine `hill_gate_index` grades its fits by. Inputs, labels and view state are not here.
 *
 * `BlockArgs` and `BlockData` are both built on this type, so a setting added here reaches the
 * workflow and the presets below together.
 */
export type FitSettings = {
  equation: Equation;
  /** The fraction of the gate range a curve must reach to count as saturated. Read by
   *  `log_fluorescence` only. */
  plateauFloor: number;
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

/** The settings only `hill_gate_index` reads, at the values the block shipped through 2.4.0. */
export const HILL_GATE_INDEX_DEFAULTS = {
  minReadsPerConcentration: 3,
  minConcentrationPoints: 5,
  r2ThresholdGood: 0.8,
  r2ThresholdFailed: 0.5,
  nMin: 0.5,
  nMax: 2.0,
  hookEffectThresholdBin: 0.2,
  hookEffectThresholdNoBin: 0.02,
  hookEffectMinReads: 20,
} as const satisfies Omit<FitSettings, "equation" | "plateauFloor">;

export const DEFAULT_EQUATION: Equation = "log_fluorescence";

/** Adams 2016's floor: two thirds of the gate range. */
export const DEFAULT_PLATEAU_FLOOR = 2 / 3;

export const PRESET_IDS = ["adams2016", "phillips2021", "starr"] as const;

export type PresetId = (typeof PRESET_IDS)[number];

export type Preset = {
  readonly id: PresetId;
  readonly label: string;
  /** The settings applying this preset writes. A setting it does not name is left as it is. */
  readonly settings: Readonly<Partial<FitSettings>>;
};

/**
 * The published analyses this block reproduces, each as the settings that reproduce it.
 *
 * A preset is not state: applying one writes the settings it names and nothing records that it
 * was the source. `matchesPreset` reads the settings back, which is what lets the drawer say a
 * preset was applied and then edited.
 */
export const PRESETS: readonly Preset[] = [
  {
    id: "adams2016",
    label: "Adams 2016",
    settings: { equation: "log_fluorescence", plateauFloor: DEFAULT_PLATEAU_FLOOR },
  },
  {
    id: "phillips2021",
    label: "Phillips 2021 (CR6261)",
    settings: { equation: "log_fluorescence", plateauFloor: 0.4 },
  },
  {
    id: "starr",
    label: "Starr 2020 / 2022",
    settings: { equation: "hill_gate_index", ...HILL_GATE_INDEX_DEFAULTS },
  },
];

export function presetById(id: PresetId): Preset {
  const preset = PRESETS.find((p) => p.id === id);
  if (preset === undefined) throw new Error(`Unknown preset: ${id}`);
  return preset;
}

/** Whether `settings` still holds every value `preset` writes. */
export function matchesPreset(preset: Preset, settings: FitSettings): boolean {
  return Object.entries(preset.settings).every(
    ([field, value]) => settings[field as keyof FitSettings] === value,
  );
}
