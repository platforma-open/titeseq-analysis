<script setup lang="ts">
import { PlAlert, PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import PageHeader from "./PageHeader.vue";

defineProps<{ title: string }>();

const app = useApp();

// Numeric-field validation surfaces inline on each PlNumberField via its
// `:error-message` prop (see useFieldValidation). The page-level alert is
// reserved for spec-based checks that have no field-level home —
// concentration column label format, sort-fraction without a bin column,
// etc. — emitted by the model's validationWarnings output.
const warnings = computed(() => app.model.outputs.validationWarnings ?? []);
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="app.model.data.defaultBlockLabel"
    :title="title"
    no-body-gutters
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
