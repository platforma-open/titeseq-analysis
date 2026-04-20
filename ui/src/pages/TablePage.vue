<script setup lang="ts">
import {
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import PageHeader from "../components/PageHeader.vue";

const app = useApp();

const warnings = computed(() => app.model.outputs.validationWarnings ?? []);

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.summaryTable,
});
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="app.model.data.defaultBlockLabel"
    title="Clonotype Fit Results"
  >
    <template #append>
      <PageHeader />
    </template>
    <PlAlert v-for="(w, i) in warnings" :key="i" :type="w.severity === 'error' ? 'error' : 'warn'">
      {{ w.message }}
    </PlAlert>
    <PlAgDataTableV2
      v-model="app.model.data.tableState"
      :settings="tableSettings"
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      no-rows-text="No clonotypes available."
    />
  </PlBlockPage>
</template>
