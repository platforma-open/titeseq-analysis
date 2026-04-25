<script setup lang="ts">
import { PlAgDataTableV2, usePlDataTableSettingsV2 } from "@platforma-sdk/ui-vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.meanBinTable,
});
</script>

<template>
  <TiteseqPage title="Mean Bin Data">
    <PlAgDataTableV2
      v-model="app.model.data.meanBinTableState"
      :settings="tableSettings"
      show-export-button
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      no-rows-text="No mean bin data available."
    />
  </TiteseqPage>
</template>
