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
  >
    <template #append>
      <PageHeader />
    </template>
    <slot />
  </PlBlockPage>
</template>
