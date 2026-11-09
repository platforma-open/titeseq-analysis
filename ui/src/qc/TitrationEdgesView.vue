<script setup lang="ts">
/**
 * How many variants the titration could not measure: bound at its lowest or highest
 * concentration. The histogram shows Kd,app by class; bounds pile up at the two edges.
 */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import { LOG_HISTOGRAM_DEFAULT_STATE } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

const app = useApp();

const graphState = computed({
  get: (): GraphMakerState =>
    app.model.data.graphStateTitrationEdges ?? LOG_HISTOGRAM_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data.graphStateTitrationEdges = state;
  },
});

const col = (name: string) => app.model.outputs.reportedPfCols?.find((p) => p.spec.name === name);

const defaultOptions = computed((): PredefinedGraphOption<"histogram">[] | undefined => {
  const kd = col("pl7.app/vdj/kd");
  const cls = col("pl7.app/vdj/affinityClass");
  if (!kd || !cls) return undefined;
  return [
    { inputName: "value", selectedSource: kd.spec },
    { inputName: "grouping", selectedSource: cls.spec },
  ];
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
