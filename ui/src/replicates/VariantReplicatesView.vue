<script setup lang="ts">
/**
 * One row per (variant, replicate): the value that replicate contributed to the mean. In
 * synonymous mode the variant is the amino-acid one, with the spread of its synonymous variants
 * inside that replicate.
 */
import { PlAgDataTableV2, usePlDataTableSettingsV2 } from "@platforma-sdk/ui-vue";
import { useApp } from "../app";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.replicateTable,
});
</script>

<template>
  <div class="titeseq-replicate-table">
    <PlAgDataTableV2
      v-model="app.model.data.replicateTableState"
      :settings="tableSettings"
      show-export-button
      not-ready-text="Configure the block in Settings on the Main page and run it."
      no-rows-text="No per-replicate values available."
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
