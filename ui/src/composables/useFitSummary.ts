import { initBlockData } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";
import { FIELDS_READ_BY, SETTINGS_FIELD_LABELS, equationLabel } from "../fitSettings";

// The block's own defaults, read from the function that seeds a fresh block.
const DEFAULTS = initBlockData();

/**
 * The line every results page carries. Fit settings sit in the collapsed tier, so a non-default
 * value is not glanceable where the numbers are read — this states which equation produced them
 * and which of the settings that equation reads differ from their defaults.
 */
export function useFitSummary() {
  const app = useApp();

  const changedSettings = computed(() => {
    const changed: string[] = [];
    for (const name of FIELDS_READ_BY[app.model.data.equation]) {
      const value = app.model.data[name];
      if (value !== DEFAULTS[name]) changed.push(`${SETTINGS_FIELD_LABELS[name]} ${value}`);
    }
    return changed;
  });

  const fitSummary = computed(() => {
    const method = `Kd,app values come from the ${equationLabel(app.model.data.equation)}.`;
    const changed = changedSettings.value;
    if (changed.length === 0) return `${method} Every setting it reads is at its default.`;
    return `${method} Changed from defaults: ${changed.join(", ")}.`;
  });

  return { fitSummary };
}
