<script setup lang="ts">
/**
 * Synonymous mode only: one row per (nt variant, replicate), the variant's own fit in that replicate,
 * beside the amino-acid variant it belongs to.
 */
import type { PlDataTableStateV2 } from "@platforma-sdk/model";
import { createPlDataTableStateV2 } from "@platforma-sdk/model";
import { PlAgDataTableV2, usePlDataTableSettingsV2 } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.replicateNtTable,
});

const tableState = computed({
  get: () => app.model.data.replicateNtTableState ?? createPlDataTableStateV2(),
  set: (value: PlDataTableStateV2) => {
    app.model.data.replicateNtTableState = value;
  },
});
</script>

<template>
  <div class="titeseq-replicate-table">
    <PlAgDataTableV2
      v-model="tableState"
      :settings="tableSettings"
      show-export-button
      not-ready-text="Configure the block in Settings on the Main page and run it."
      no-rows-text="No per-replicate nucleotide values available."
    />
  </div>
</template>

<style scoped>
.titeseq-replicate-table {
  flex: 1;
  display: flex;
  flex-direction: column;
  padding: 8px 16px 16px;
}
</style>
