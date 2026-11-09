<script setup lang="ts">
/** Reads behind each fit, one box per class: whether fits fail for lack of reads. */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import { READ_DEPTH_DEFAULT_STATE } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

const app = useApp();

const graphState = computed({
  get: (): GraphMakerState => app.model.data.graphStateReadsByClass ?? READ_DEPTH_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data.graphStateReadsByClass = state;
  },
});

const col = (name: string) => app.model.outputs.reportedPfCols?.find((p) => p.spec.name === name);

const defaultOptions = computed((): PredefinedGraphOption<"discrete">[] | undefined => {
  const reads = col("pl7.app/vdj/fitReads");
  const cls = col("pl7.app/vdj/affinityClass");
  if (!reads || !cls) return undefined;
  return [
    { inputName: "y", selectedSource: reads.spec },
    { inputName: "primaryGrouping", selectedSource: cls.spec },
  ];
});
</script>

<template>
  <GraphMaker
    v-model="graphState"
    chart-type="discrete"
    :data-state-key="app.model.outputs.reportedPfHandle"
    :p-frame="app.model.outputs.reportedPfHandle"
    :default-options="defaultOptions"
    :status-text="{
      noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
    }"
  />
</template>
