<script setup lang="ts">
/**
 * Synonymous mode only: one row per nucleotide variant, its own fit beside its class and parent
 * delta. `ntVariantTable` is an `outputWithStatus`, so the table owns its empty states.
 */
import { PlAgDataTableV2, usePlDataTableSettingsV2 } from "@platforma-sdk/ui-vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.ntVariantTable,
});
</script>

<template>
  <TiteseqPage title="Nucleotide Variant Fit Results" mode="table">
    <PlAgDataTableV2
      v-model="app.model.data.ntTableState"
      :settings="tableSettings"
      show-export-button
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      no-rows-text="No nucleotide variants available."
    />
  </TiteseqPage>
</template>
