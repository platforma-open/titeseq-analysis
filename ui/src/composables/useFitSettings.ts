import type {
  Equation,
  FitSettings,
  PresetId,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  PRESETS,
  matchesPreset,
  presetById,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

/**
 * Writing fit settings into block data, and preset selection over the model's one preset table.
 *
 * A preset is a starting point: applying it writes the real fields, and `appliedPreset` records
 * only which preset was the source. Whether the settings still match is read back from the
 * settings themselves, so a later edit shows as "(modified)" instead of being reverted.
 */
export function useFitSettings() {
  const app = useApp();

  /**
   * Writes the equation and the Affinity vs Fit Quality y-axis scale. The y column is the SE under
   * log_fluorescence (log scale) and hillPlotPosition under hill_gate_index (linear: it holds -1).
   */
  function setEquation(equation: Equation): void {
    app.model.data.equation = equation;
    const graph = app.model.data.graphStateAffinityVsFit;
    graph.axesSettings = {
      ...graph.axesSettings,
      axisY: {
        ...graph.axesSettings?.axisY,
        scale: equation === "log_fluorescence" ? "log" : "linear",
      },
    };
  }

  /** Writes every setting the argument names, and leaves the rest as they are. */
  function applySettings({ equation, ...rest }: Readonly<Partial<FitSettings>>): void {
    Object.assign(app.model.data, rest);
    if (equation !== undefined) setEquation(equation);
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

  return { presetOptions, selectedPreset, isModified, setEquation };
}
