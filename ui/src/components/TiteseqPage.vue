<script setup lang="ts">
import { PlBlockPage } from "@platforma-sdk/ui-vue";
import { useApp } from "../app";
import PageHeader from "./PageHeader.vue";

defineProps<{ title: string }>();

const app = useApp();
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
    <slot />
  </PlBlockPage>
</template>

<style scoped>
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
