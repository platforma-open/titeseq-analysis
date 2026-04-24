<script setup lang="ts">
import { PlAlert, PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import { useFieldValidation } from "../composables/useFieldValidation";
import PageHeader from "./PageHeader.vue";

defineProps<{ title: string }>();

const app = useApp();

// Local field-level validation runs synchronously against app.model.data, so
// the page-level alert updates the moment the user types — without waiting on
// the mutation→server→patch round-trip that drives outputs.validationWarnings.
const { warnings: localWarnings } = useFieldValidation();

// Server-side validationWarnings still covers spec-based checks (concentration
// column label, antigen/target pairing) that need ctx.resultPool. Merge both
// and dedupe by message so the user sees the full picture without duplicates.
const warnings = computed(() => {
  const serverWarnings = app.model.outputs.validationWarnings ?? [];
  const seen = new Set<string>();
  return [...localWarnings.value, ...serverWarnings].filter((w) => {
    if (seen.has(w.message)) return false;
    seen.add(w.message);
    return true;
  });
});
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
