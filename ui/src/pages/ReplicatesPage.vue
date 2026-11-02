<script setup lang="ts">
/**
 * Replicate mode only: each replicate's own Kd,app for each variant, as a table and as a chart
 * that puts the replicates of one variant side by side. The combined value lives in the main
 * table; this page is where the spread behind it is read.
 *
 * `replicateTable` is an `outputWithStatus`, so the table owns its empty states.
 */
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { PlAgDataTableV2, PlBtnGroup, usePlDataTableSettingsV2 } from "@platforma-sdk/ui-vue";
import { computed, ref } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const VIEWS = [
  { value: "table" as const, label: "Table" },
  { value: "chart" as const, label: "Chart" },
];
const view = ref<"table" | "chart">("table");

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.replicateTable,
});

const REPLICATE_AXIS = "pl7.app/vdj/replicate";

const defaultOptions = computed((): PredefinedGraphOption<"scatterplot">[] | undefined => {
  const pCols = app.model.outputs.replicatePfCols;
  if (!pCols) return undefined;

  const kd = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/kd");
  if (!kd) return undefined;

  const replicateAxis = kd.spec.axesSpec.find((a) => a.name === REPLICATE_AXIS);
  // The variant axis is whichever one is not the replicate: clonotypeKey, scClonotypeKey and
  // variantKey all resolve here without naming each of them.
  const variantAxis = kd.spec.axesSpec.find((a) => a.name !== REPLICATE_AXIS);
  if (!replicateAxis || !variantAxis) return undefined;

  const options: PredefinedGraphOption<"scatterplot">[] = [
    { inputName: "x", selectedSource: replicateAxis },
    { inputName: "y", selectedSource: kd.spec },
    { inputName: "facetBy", selectedSource: variantAxis },
    { inputName: "grouping", selectedSource: replicateAxis },
  ];

  const state = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/kdReplicateState");
  if (state) {
    options.push({
      inputName: "filters",
      selectedSource: state.spec,
      filterType: "equals",
      selectedFilterValues: ["measured"],
    });
  }
  return options;
});
</script>

<template>
  <TiteseqPage title="Replicates" :mode="view === 'table' ? 'table' : 'graph'">
    <PlBtnGroup v-model="view" :options="VIEWS" compact class="titeseq-replicate-view" />
    <PlAgDataTableV2
      v-if="view === 'table'"
      v-model="app.model.data.replicateTableState"
      :settings="tableSettings"
      show-export-button
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      no-rows-text="No per-replicate values available."
    />
    <GraphMaker
      v-else
      v-model="app.model.data.graphStateReplicates"
      chart-type="scatterplot"
      :data-state-key="app.model.outputs.replicatePfHandle"
      :p-frame="app.model.outputs.replicatePfHandle"
      :default-options="defaultOptions"
      default-palette="bright"
      :status-text="{
        noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
      }"
    />
  </TiteseqPage>
</template>

<style scoped>
.titeseq-replicate-view {
  margin: 8px 16px;
  max-width: 240px;
}
</style>
