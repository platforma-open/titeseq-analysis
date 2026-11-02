<script setup lang="ts">
import { PlAgDataTableV2, usePlDataTableSettingsV2 } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";
import { useNoun } from "../composables/useNoun";

const app = useApp();
const { noun, Noun } = useNoun();

// Synonymous mode shows one row per amino-acid variant.
const title = computed(() =>
  app.model.outputs.synonymous === true
    ? "Amino Acid Variant Fit Results"
    : `${Noun.value} Fit Results`,
);
const noRowsText = computed(() =>
  app.model.outputs.synonymous === true
    ? "No amino acid variants available."
    : `No ${noun.value}s available.`,
);

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.summaryTable,
});
</script>

<template>
  <TiteseqPage :title="title" mode="table">
    <PlAgDataTableV2
      v-model="app.model.data.tableState"
      :settings="tableSettings"
      show-export-button
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      :no-rows-text="noRowsText"
    />
  </TiteseqPage>
</template>
