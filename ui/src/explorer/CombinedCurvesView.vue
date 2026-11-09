<script setup lang="ts">
/**
 * The QC page's titration curves with the replicates combined: one panel per variant, its observed
 * signal and its fitted curve at every concentration, each a point at the mean over replicates with
 * a standard-error bar. No lines: a category chart draws the same layers for both series. A category chart, the only chart that draws error bars, so
 * the concentrations sit evenly spaced. It takes one value column, so the two come stacked under a
 * series axis. The multi-select picks the variants by their affinity class across the replicates.
 */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import {
  EXPLORER_CURVE_CLASSES_DEFAULT,
  EXPLORER_CURVES_DEFAULT_STATE,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import { PlDropdownMulti } from "@platforma-sdk/ui-vue";
// The package root exports only `GraphMaker` and types; this runtime function is reachable only
// through the `./*` wildcard, which maps literally and so needs the `.js`.
import { createDefaultMapping } from "@milaboratories/graph-maker/dist/dataBindAes.js";
import { computed, watch } from "vue";
import { useApp } from "../app";
import { axisId } from "../composables/useConcentrationOrder";
import { resetExplorerCurves } from "./resetExplorerCurves";

const app = useApp();

const CONCENTRATION_AXIS = "pl7.app/vdj/concentration";
const REPLICATE_AXIS = "pl7.app/vdj/replicate";
const SERIES_AXIS = "pl7.app/vdj/signalSeries";
const SERIES_ORDER = ["Observed", "Fitted"];

const cols = computed(() => app.model.outputs.titrationCurvesPfCols);
const signal = computed(() =>
  cols.value?.find((p) => p.spec.name === "pl7.app/vdj/titrationSignal"),
);
const seriesAxis = computed(() => signal.value?.spec.axesSpec.find((a) => a.name === SERIES_AXIS));
// The variant's class across the replicates: the summary column, on the variant axis alone.
const affinityClass = computed(() =>
  cols.value?.find(
    (p) => p.spec.name === "pl7.app/vdj/affinityClass" && p.spec.axesSpec.length === 1,
  ),
);
// The concentration as a label that sorts as text in numeric order: a faceted category chart
// sorts its groups as text and ignores a stored order.
const concentrationLabel = computed(() =>
  cols.value?.find((p) => p.spec.name === "pl7.app/vdj/concentrationLabel"),
);

/** The classes the equation can rate, in its declared order. */
const classOptions = computed(() => {
  const declared = affinityClass.value?.spec.annotations?.["pl7.app/discreteValues"];
  const values: unknown = declared ? JSON.parse(declared) : [];
  return Array.isArray(values) ? values.map((v) => ({ label: String(v), value: String(v) })) : [];
});

/**
 * The state must exist in data before the binding reads it: a getter that returns a new default
 * object per read never lets the chart settle, and it renders blank.
 */
app.model.data.graphStateExplorerCurves ??= { ...EXPLORER_CURVES_DEFAULT_STATE };

const storedState = computed({
  get: () => app.model.data.graphStateExplorerCurves,
  set: (state: GraphMakerState | undefined) => {
    app.model.data.graphStateExplorerCurves = state;
  },
});

/**
 * Observed first, fitted dotted where a line layer is turned on. The chart owns the series entry: it deletes it while the frame has
 * no values and rebuilds it with its own defaults. So the entry is corrected whenever it differs,
 * never created ahead of it.
 */
watch(
  [
    seriesAxis,
    () => seriesAxis.value && storedState.value?.dataBindAes?.[axisId(seriesAxis.value)],
  ],
  ([axis, entry]) => {
    const current = storedState.value;
    if (!axis || !current || entry?.type !== "categorical") return;
    const fittedDotted = entry.mapping.Fitted?.aes.lineShape === "dotted";
    const inOrder = SERIES_ORDER.every((v, i) => entry.order[i] === v);
    if (fittedDotted && inOrder) return;
    const mapping = createDefaultMapping([...SERIES_ORDER], "bright");
    mapping.mapping.Observed.aes.lineShape = "solid";
    mapping.mapping.Fitted.aes.lineShape = "dotted";
    storedState.value = {
      ...current,
      dataBindAes: { ...current.dataBindAes, [axisId(axis)]: mapping },
    };
  },
  { immediate: true },
);

const graphState = computed({
  get: (): GraphMakerState => storedState.value ?? EXPLORER_CURVES_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    storedState.value = state;
  },
});

/**
 * The picked classes reach the chart as its default filter. The chart applies its defaults once,
 * so a new pick starts the chart again from the default state, and the key remounts it.
 */
const classes = computed({
  get: () => app.model.data.explorerCurveClasses ?? EXPLORER_CURVE_CLASSES_DEFAULT,
  set: (picked: string[]) => {
    app.model.data.explorerCurveClasses = picked;
    resetExplorerCurves(app.model.data);
  },
});
const chartKey = computed(() => classes.value.join("|"));

const defaultOptions = computed((): PredefinedGraphOption<"discrete">[] | undefined => {
  const y = signal.value;
  const concentration = concentrationLabel.value?.spec;
  const variantAxis = y?.spec.axesSpec.find(
    (a) => a.name !== CONCENTRATION_AXIS && a.name !== REPLICATE_AXIS && a.name !== SERIES_AXIS,
  );
  const series = seriesAxis.value;
  if (!y || !concentration || !variantAxis || !series || !affinityClass.value) return undefined;
  return [
    { inputName: "y", selectedSource: y.spec },
    { inputName: "primaryGrouping", selectedSource: concentration },
    { inputName: "secondaryGrouping", selectedSource: series },
    { inputName: "facetBy", selectedSource: variantAxis },
    {
      inputName: "filters",
      selectedSource: affinityClass.value.spec,
      filterType: "equals",
      selectedFilterValues: classes.value,
    },
  ];
});
</script>

<template>
  <div class="titeseq-combined-curves">
    <PlDropdownMulti
      v-model="classes"
      :options="classOptions"
      label="Affinity class (combined over replicates)"
      class="titeseq-combined-curves-classes"
    />
    <GraphMaker
      :key="chartKey"
      v-model="graphState"
      chart-type="discrete"
      :p-frame="app.model.outputs.titrationCurvesPf"
      :default-options="defaultOptions"
      default-palette="bright"
      :status-text="{
        noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
      }"
    />
  </div>
</template>

<style scoped>
.titeseq-combined-curves {
  flex: 1;
  display: flex;
  flex-direction: column;
  min-height: 480px;
}
.titeseq-combined-curves-classes {
  max-width: 420px;
  margin: 8px 16px;
}
</style>
