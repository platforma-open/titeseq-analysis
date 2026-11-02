import type { Equation } from "@platforma-open/platforma-open.titeseq-analysis.kind";

/**
 * Every setting that shapes the fit — the two knobs the equations are chosen and tuned with, and
 * the nine `hill_gate_index` grades its fits by. Inputs, labels and view state are not here.
 *
 * The `hill_gate_index` branch of `BlockArgs` is built on this type, and `BlockData` on
 * `FitSettingsDraft` below, so a setting added here reaches the drawer, the args and the presets
 * together. Every value is present.
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
export type LogFluorescenceSetting = "plateauFloor" | "r2GateEnabled" | "r2GateMin";

/**
 * The settings as the block stores them, which is what the drawer edits: a number field the
 * scientist cleared reads `undefined`. `BlockData` is built on this type, and `args()` is the
 * single boundary that refuses an incomplete set and so disables Run.
 *
 * `equation` and `r2GateEnabled` stay required: a dropdown and a checkbox have no clear button.
 */
export type FitSettingsDraft = Partial<Omit<FitSettings, "equation" | "r2GateEnabled">> &
  Pick<FitSettings, "equation" | "r2GateEnabled">;

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
} as const satisfies Omit<FitSettings, "equation" | LogFluorescenceSetting>;

export const DEFAULT_EQUATION: Equation = "log_fluorescence";

/** Adams 2016's floor: two thirds of the gate range. */
export const DEFAULT_PLATEAU_FLOOR = 2 / 3;

/** Off: no rule unverified against a paper's own code sits in the default path. The threshold is
 *  Phillips 2021's, which its preset turns on. */
export const DEFAULT_R2_GATE = { r2GateEnabled: false, r2GateMin: 0.8 } as const;

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
    settings: {
      equation: "log_fluorescence",
      plateauFloor: DEFAULT_PLATEAU_FLOOR,
      r2GateEnabled: false,
    },
  },
  {
    id: "phillips2021",
    label: "Phillips 2021 (CR6261)",
    // The authors drop an in-range fit with R² below 0.8 (klawrence26/bnab-landscapes).
    settings: {
      equation: "log_fluorescence",
      plateauFloor: 0.4,
      r2GateEnabled: true,
      r2GateMin: 0.8,
    },
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
export function matchesPreset(preset: Preset, settings: FitSettingsDraft): boolean {
  return Object.entries(preset.settings).every(
    ([field, value]) => settings[field as keyof FitSettings] === value,
  );
}
