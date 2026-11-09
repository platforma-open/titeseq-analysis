<script setup lang="ts">
/**
 * The titration curves with the replicates combined, on a log concentration axis: one panel per
 * variant, the observed mean over replicates as points and the replicates' mean curve as a line.
 * Both average only the replicates whose fit drew a curve. No error bars: a scatterplot has none,
 * though the standard error is in the frame (`pl7.app/vdj/combinedMeanBinSe`).
 */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import {
  EXPLORER_CURVE_CLASSES_DEFAULT,
  EXPLORER_SCATTER_DEFAULT_STATE,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import { PlDropdownMulti } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import { resetExplorerCurves } from "./resetExplorerCurves";

const app = useApp();

const CONCENTRATION_AXIS = "pl7.app/vdj/concentration";

const cols = computed(() => app.model.outputs.titrationCurvesPfCols);
const byName = (name: string) => cols.value?.find((p) => p.spec.name === name);
// The variant's class across the replicates: the summary column, on the variant axis alone.
const affinityClass = computed(() =>
  cols.value?.find(
    (p) => p.spec.name === "pl7.app/vdj/affinityClass" && p.spec.axesSpec.length === 1,
  ),
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
app.model.data.graphStateExplorerScatter ??= { ...EXPLORER_SCATTER_DEFAULT_STATE };

const graphState = computed({
  get: (): GraphMakerState =>
    app.model.data.graphStateExplorerScatter ?? EXPLORER_SCATTER_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data.graphStateExplorerScatter = state;
  },
});

/** Shared with the mean ± SE tab. A new pick starts both charts again from their defaults. */
const classes = computed({
  get: () => app.model.data.explorerCurveClasses ?? EXPLORER_CURVE_CLASSES_DEFAULT,
  set: (picked: string[]) => {
    app.model.data.explorerCurveClasses = picked;
    resetExplorerCurves(app.model.data);
  },
});
const chartKey = computed(() => classes.value.join("|"));

const defaultOptions = computed((): PredefinedGraphOption<"scatterplot">[] | undefined => {
  const mean = byName("pl7.app/vdj/combinedMeanBin");
  const curve = byName("pl7.app/vdj/combinedFittedMeanBin");
  const concValue = byName("pl7.app/vdj/concentrationValue");
  const variantAxis = mean?.spec.axesSpec.find((a) => a.name !== CONCENTRATION_AXIS);
  if (!mean || !curve || !concValue || !variantAxis || !affinityClass.value) return undefined;
  return [
    { inputName: "x", selectedSource: concValue.spec },
    { inputName: "y", selectedSource: mean.spec },
    { inputName: "facetBy", selectedSource: variantAxis },
    { inputName: "additionalCurves", selectedSource: curve.spec },
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
  <div class="titeseq-scatter-curves">
    <PlDropdownMulti
      v-model="classes"
      :options="classOptions"
      label="Affinity class (combined over replicates)"
      class="titeseq-scatter-curves-classes"
    />
    <GraphMaker
      :key="chartKey"
      v-model="graphState"
      chart-type="scatterplot"
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
.titeseq-scatter-curves {
  flex: 1;
  display: flex;
  flex-direction: column;
  min-height: 480px;
}
.titeseq-scatter-curves-classes {
  max-width: 420px;
  margin: 8px 16px;
}
</style>
