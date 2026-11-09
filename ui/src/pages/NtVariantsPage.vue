<script setup lang="ts">
/**
 * Synonymous mode only: one row per nucleotide variant, its own fit beside its class and parent
 * delta. `ntVariantTable` is an `outputWithStatus`, so the table owns its empty states.
 */
import { PlAgDataTableV2, PlTabs, usePlDataTableSettingsV2 } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.ntVariantTable,
});

// The model filters the rows; the synonymous variants are the ones the synonymous graphs list.
type Shown = "all" | "synonymous";
const shown = computed({
  get: (): Shown => (app.model.data.ntSynonymousOnly === true ? "synonymous" : "all"),
  set: (value: Shown) => {
    app.model.data.ntSynonymousOnly = value === "synonymous";
  },
});
const shownOptions: { label: string; value: Shown }[] = [
  { label: "All nucleotide variants", value: "all" },
  { label: "Synonymous variants only", value: "synonymous" },
];
</script>

<template>
  <TiteseqPage title="Nucleotide Variant Fit Results" mode="table">
    <PlTabs v-model="shown" :options="shownOptions" class="titeseq-nt-shown" />
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
  margin-bottom: 8px;
}
</style>
