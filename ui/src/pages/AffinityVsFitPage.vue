<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import PageHeader from "../components/PageHeader.vue";

const app = useApp();

const defaultOptions = computed((): PredefinedGraphOption<"scatterplot">[] | undefined => {
  const pCols = app.model.outputs.summaryPfCols;
  if (!pCols) return undefined;

  const kdPlot = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/kdPlotPosition");
  const hillPlot = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/hillPlotPosition",
  );
  const affinityClass = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/affinityClass",
  );
  const failureReason = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/fitFailureReason",
  );
  if (!kdPlot || !hillPlot || !affinityClass) return undefined;

  const options: PredefinedGraphOption<"scatterplot">[] = [
    { inputName: "x", selectedSource: kdPlot.spec },
    { inputName: "y", selectedSource: hillPlot.spec },
    { inputName: "grouping", selectedSource: affinityClass.spec },
  ];
  if (failureReason) {
    options.push({ inputName: "shape", selectedSource: failureReason.spec });
  }
  return options;
});
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.args.customBlockLabel"
    :subtitle-placeholder="app.model.args.defaultBlockLabel"
    title="Affinity vs Fit Quality"
  >
    <template #append>
      <PageHeader />
    </template>
    <GraphMaker
      v-model="app.model.ui.graphStateAffinityVsFit"
      chart-type="scatterplot"
      :data-state-key="app.model.outputs.summaryPfHandle"
      :p-frame="app.model.outputs.summaryPfHandle"
      :default-options="defaultOptions"
      :status-text="{ noPframe: { title: 'Configure inputs in Settings and run the block.' } }"
    />
  </PlBlockPage>
</template>
