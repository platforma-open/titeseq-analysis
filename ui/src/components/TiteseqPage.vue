<script setup lang="ts">
import { PlAlert, PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import PageHeader from "./PageHeader.vue";

defineProps<{ title: string; mode: "graph" | "table" }>();

const app = useApp();

// Page-level alerts cover spec-based checks with no field-level home, such as a
// sort fraction without a bin column. Issues scoped to "settings" show in the
// settings drawer only. Numeric-field violations show inline on each
// PlNumberField (see useFieldValidation).
const warnings = computed(() =>
  (app.model.outputs.validationWarnings ?? []).filter((w) => w.scope === "everywhere"),
);
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="app.model.data.defaultBlockLabel"
    :title="title"
    :no-body-gutters="mode === 'graph'"
  >
    <template #append>
      <PageHeader />
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
