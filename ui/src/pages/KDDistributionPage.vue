<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import PageHeader from "../components/PageHeader.vue";

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
      selectedFilterValues: ["Good", "Partial"],
    },
  ];
});
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="app.model.data.defaultBlockLabel"
    title="K_D Distribution"
  >
    <template #append>
      <PageHeader />
    </template>
    <GraphMaker
      v-model="app.model.data.graphStateKDHistogram"
      chart-type="histogram"
      :data-state-key="app.model.outputs.summaryPfHandle"
      :p-frame="app.model.outputs.summaryPfHandle"
      :default-options="defaultOptions"
      :status-text="{
        noPframe: { title: 'Configure inputs on the Overview tab and run the block.' },
      }"
    />
  </PlBlockPage>
</template>
