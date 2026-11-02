<script setup lang="ts">
import { PlAlert, PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import PageHeader from "./PageHeader.vue";

// Only the main page offers Settings.
defineProps<{ title: string; mode: "graph" | "table"; showSettings?: boolean }>();

const app = useApp();

// Checks with no field of their own. Numeric-field errors show on each PlNumberField.
const warnings = computed(() => app.model.outputs.validationWarnings ?? []);
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="app.model.data.defaultBlockLabel"
    :title="title"
    :no-body-gutters="mode === 'graph'"
  >
    <template #append>
      <PageHeader :show-settings="showSettings" />
    </template>
    <PlAlert
      v-for="(w, i) in warnings"
      :key="i"
      :type="w.severity === 'error' ? 'error' : 'warn'"
      class="titeseq-page-alert"
    >
      {{ w.message }}
    </PlAlert>
    <slot />
  </PlBlockPage>
</template>

<style scoped>
.titeseq-page-alert {
  margin: 8px 16px 0;
}
:deep(.graph-maker .chart_title),
:deep(.graph-maker .chart_titleEdit) {
  display: none;
}
:deep(.graph-maker .chart_header) {
  margin-bottom: 0;
}
:deep(.graph-maker .chart_container) {
  padding-top: 0;
}
</style>
