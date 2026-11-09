<script setup lang="ts">
/**
 * Synonymous mode only: one row per nucleotide variant, its own fit beside its class and parent
 * delta. `ntVariantTable` is an `outputWithStatus`, so the table owns its empty states.
 */
import { PlAgDataTableV2, PlBtnGroup, usePlDataTableSettingsV2 } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.ntVariantTable,
});

// The model filters the rows; the synonymous copies are the ones the synonymous graphs list.
const synonymousOnly = computed({
  get: () => app.model.data.ntSynonymousOnly ?? false,
  set: (value: boolean) => {
    app.model.data.ntSynonymousOnly = value;
  },
});
const shownOptions = [
  { label: "All nucleotide variants", value: false },
  { label: "Synonymous copies only", value: true },
];
</script>

<template>
  <TiteseqPage title="Nucleotide Variant Fit Results" mode="table">
    <PlBtnGroup v-model="synonymousOnly" :options="shownOptions" class="titeseq-nt-shown" />
    <PlAgDataTableV2
      v-model="app.model.data.ntTableState"
      :settings="tableSettings"
      show-export-button
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      no-rows-text="No nucleotide variants available."
    />
  </TiteseqPage>
</template>

<style scoped>
.titeseq-nt-shown {
  align-self: flex-start;
  margin-bottom: 8px;
}
</style>
