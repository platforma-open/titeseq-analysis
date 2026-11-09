<script setup lang="ts">
/** Parent mode only: a histogram of each variant's log10 Kd,app vs Parent. */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import {
  PARENT_DELTA_COL,
  VS_PARENT_DEFAULT_STATE,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

const app = useApp();

const graphState = computed({
  get: (): GraphMakerState => app.model.data.graphStateVsParent ?? VS_PARENT_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data.graphStateVsParent = state;
  },
});

const col = (name: string) => app.model.outputs.reportedPfCols?.find((p) => p.spec.name === name);

const defaultOptions = computed((): PredefinedGraphOption<"histogram">[] | undefined => {
  const delta = col(PARENT_DELTA_COL);
  if (!delta) return undefined;
  return [{ inputName: "value", selectedSource: delta.spec }];
});
</script>

<template>
  <GraphMaker
    v-model="graphState"
    chart-type="histogram"
    :data-state-key="app.model.outputs.reportedPfHandle"
    :p-frame="app.model.outputs.reportedPfHandle"
    :default-options="defaultOptions"
    :status-text="{
      noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
    }"
  />
</template>
