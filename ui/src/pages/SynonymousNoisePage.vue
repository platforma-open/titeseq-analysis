<script setup lang="ts">
/**
 * Synonymous mode only: the spread of log10 Kd,app across synonymous copies of one amino-acid
 * variant. The copies encode the same protein, so the spread is the measurement noise.
 */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import { NOISE_HISTOGRAM_DEFAULT_STATE } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const graphState = computed({
  get: (): GraphMakerState =>
    app.model.data.graphStateSynonymousNoise ?? NOISE_HISTOGRAM_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data.graphStateSynonymousNoise = state;
  },
});

const defaultOptions = computed((): PredefinedGraphOption<"histogram">[] | undefined => {
  const spread = app.model.outputs.reportedPfCols?.find(
    (p) => p.spec.name === "pl7.app/vdj/kdLog10SynonymSd",
  );
  if (!spread) return undefined;
  return [{ inputName: "value", selectedSource: spread.spec }];
});
</script>

<template>
  <TiteseqPage title="Synonymous Noise" mode="graph">
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
  </TiteseqPage>
</template>
