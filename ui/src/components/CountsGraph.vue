<script setup lang="ts">
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";

// A bar chart of one run-quality count column: one bar per category, or with replicates one bar
// per replicate stacked by category.
const props = defineProps<{ columnName: string }>();
const state = defineModel<GraphMakerState | undefined>({ required: true });

const app = useApp();

const REPLICATE_AXIS = "pl7.app/vdj/replicate";

const column = computed(() =>
  app.model.outputs.runQualityPfCols?.find(
    (p: PColumnIdAndSpec) => p.spec.name === props.columnName,
  ),
);
const replicateAxis = computed(() =>
  column.value?.spec.axesSpec.find((a) => a.name === REPLICATE_AXIS),
);

const BARS: GraphMakerState = {
  title: "",
  template: "bar",
  currentTab: null,
  layersSettings: { bar: {} },
};
const STACKED_BARS: GraphMakerState = {
  title: "",
  template: "stackedBar",
  currentTab: null,
  layersSettings: { stackedBar: {} },
};

const graphState = computed<GraphMakerState>({
  get: () => {
    if (state.value !== undefined) return state.value;
    return replicateAxis.value ? STACKED_BARS : BARS;
  },
  set: (value: GraphMakerState) => {
    state.value = value;
  },
});

const defaultOptions = computed((): PredefinedGraphOption<"discrete">[] | undefined => {
  const spec = column.value?.spec;
  if (!spec) return undefined;
  const category = spec.axesSpec.find((a) => a.name !== REPLICATE_AXIS);
  if (!category) return undefined;
  const replicate = replicateAxis.value;
  return [
    { inputName: "y", selectedSource: spec },
    ...(replicate
      ? [
          { inputName: "primaryGrouping" as const, selectedSource: replicate },
          { inputName: "secondaryGrouping" as const, selectedSource: category },
        ]
      : [{ inputName: "primaryGrouping" as const, selectedSource: category }]),
  ];
});
</script>

<template>
  <GraphMaker
    v-if="column"
    v-model="graphState"
    chart-type="discrete"
    :data-state-key="app.model.outputs.runQualityPf"
    :p-frame="app.model.outputs.runQualityPf"
    :default-options="defaultOptions"
    :status-text="{
      noPframe: { title: 'Configure the block in Settings on the Main page and run it.' },
    }"
  />
</template>
