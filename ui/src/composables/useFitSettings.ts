import type { FitSettings, PresetId } from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  PRESETS,
  matchesPreset,
  presetById,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

type PartialSettings = Readonly<Partial<FitSettings>>;

/**
 * Writing fit settings into block data, and preset selection over the model's one preset table.
 *
 * A preset is a starting point: applying it writes the real fields, and `appliedPreset` records
 * only which preset was the source. Whether the settings still match is read back from the
 * settings themselves, so a later edit shows as "(modified)" instead of being reverted.
 */
export function useFitSettings() {
  const app = useApp();

  // Exhaustive over FitSettings: a setting added to the model fails to compile until it can be
  // written here. Each writer skips the setting its source leaves out.
  const writers: Record<keyof FitSettings, (settings: PartialSettings) => void> = {
    equation: ({ equation }) => {
      if (equation !== undefined) app.model.data.equation = equation;
    },
    plateauFloor: ({ plateauFloor }) => {
      if (plateauFloor !== undefined) app.model.data.plateauFloor = plateauFloor;
    },
    minReadsPerConcentration: ({ minReadsPerConcentration }) => {
      if (minReadsPerConcentration !== undefined)
        app.model.data.minReadsPerConcentration = minReadsPerConcentration;
    },
    minConcentrationPoints: ({ minConcentrationPoints }) => {
      if (minConcentrationPoints !== undefined)
        app.model.data.minConcentrationPoints = minConcentrationPoints;
    },
    r2ThresholdGood: ({ r2ThresholdGood }) => {
      if (r2ThresholdGood !== undefined) app.model.data.r2ThresholdGood = r2ThresholdGood;
    },
    r2ThresholdFailed: ({ r2ThresholdFailed }) => {
      if (r2ThresholdFailed !== undefined) app.model.data.r2ThresholdFailed = r2ThresholdFailed;
    },
    nMin: ({ nMin }) => {
      if (nMin !== undefined) app.model.data.nMin = nMin;
    },
    nMax: ({ nMax }) => {
      if (nMax !== undefined) app.model.data.nMax = nMax;
    },
    hookEffectThresholdBin: ({ hookEffectThresholdBin }) => {
      if (hookEffectThresholdBin !== undefined)
        app.model.data.hookEffectThresholdBin = hookEffectThresholdBin;
    },
    hookEffectThresholdNoBin: ({ hookEffectThresholdNoBin }) => {
      if (hookEffectThresholdNoBin !== undefined)
        app.model.data.hookEffectThresholdNoBin = hookEffectThresholdNoBin;
    },
    hookEffectMinReads: ({ hookEffectMinReads }) => {
      if (hookEffectMinReads !== undefined) app.model.data.hookEffectMinReads = hookEffectMinReads;
    },
  };

  /** Writes every setting the argument names, and leaves the rest as they are. */
  function applySettings(settings: PartialSettings): void {
    for (const write of Object.values(writers)) write(settings);
  }

  const isModified = computed(() => {
    const id = app.model.data.appliedPreset;
    if (id === undefined) return false;
    return !matchesPreset(presetById(id), app.model.data);
  });

  const presetOptions = computed(() =>
    PRESETS.map((preset) => ({
      value: preset.id,
      label:
        preset.id === app.model.data.appliedPreset && isModified.value
          ? `${preset.label} (modified)`
          : preset.label,
    })),
  );

  const selectedPreset = computed<PresetId | undefined>({
    get: () => app.model.data.appliedPreset,
    set: (id) => {
      if (id !== undefined) applySettings(presetById(id).settings);
      app.model.data.appliedPreset = id;
    },
  });

  return { applySettings, presetOptions, selectedPreset, isModified };
}
