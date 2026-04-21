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
  const hillCoef = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/hillCoefficient",
  );
  const affinityClass = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/affinityClass",
  );
  const failureReason = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/fitFailureReason",
  );
  if (!kdPlot || !hillPlot || !failureReason) return undefined;

  const options: PredefinedGraphOption<"scatterplot">[] = [
    { inputName: "x", selectedSource: kdPlot.spec },
    { inputName: "y", selectedSource: hillPlot.spec },
    { inputName: "grouping", selectedSource: failureReason.spec },
  ];
  if (hillCoef) {
    options.push({
      inputName: "filters",
      selectedSource: hillCoef.spec,
      filterType: "range",
      selectedFilterRange: { min: 0 },
    });
  }
  if (affinityClass) {
    options.push({
      inputName: "filters",
      selectedSource: affinityClass.spec,
      filterType: "equals",
      selectedFilterValues: ["Failed"],
    });
  }
  return options;
});
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="app.model.data.defaultBlockLabel"
    title="Affinity vs Fit Quality"
  >
    <template #append>
      <PageHeader />
    </template>
    <GraphMaker
      v-model="app.model.data.graphStateAffinityVsFit"
      chart-type="scatterplot"
      :data-state-key="app.model.outputs.summaryPfHandle"
      :p-frame="app.model.outputs.summaryPfHandle"
      :default-options="defaultOptions"
      :status-text="{
        noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
      }"
    />
  </PlBlockPage>
</template>
