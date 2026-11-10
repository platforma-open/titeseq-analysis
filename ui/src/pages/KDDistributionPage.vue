<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const defaultOptions = computed((): PredefinedGraphOption<"histogram">[] | undefined => {
  const pCols = app.model.outputs.summaryPfCols;
  if (!pCols) return undefined;

  const kd = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/kd");
  const affinityClass = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/affinityClass",
  );
  if (!kd || !affinityClass) return undefined;

  return [
    { inputName: "value", selectedSource: kd.spec },
    {
      inputName: "filters",
      selectedSource: affinityClass.spec,
      selectedFilterValues: ["Good"],
    },
  ];
});
</script>

<template>
  <TiteseqPage title="Kd Distribution" mode="graph">
    <GraphMaker
      v-model="app.model.data.graphStateKDHistogram"
      chart-type="histogram"
      :data-state-key="app.model.outputs.summaryPfHandle"
      :p-frame="app.model.outputs.summaryPfHandle"
      :default-options="defaultOptions"
      :status-text="{
        noPframe: { title: 'Configure the block in Settings on the Main page and run it.' },
      }"
    />
  </TiteseqPage>
</template>
