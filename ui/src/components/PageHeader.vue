<script setup lang="ts">
import { plRefsEqual } from "@platforma-sdk/model";
import { PlBtnGhost, PlLogView, PlMaskIcon24, PlSlideModal } from "@platforma-sdk/ui-vue";
import { onMounted, ref, watchEffect } from "vue";
import { useApp } from "../app";
import SettingsDrawer from "./SettingsDrawer.vue";

const app = useApp();
const logOpen = ref(false);

onMounted(() => {
  if (app.model.data.abundanceRef === undefined && !app.model.data.settingsOpen) {
    app.model.data.settingsOpen = true;
  }
});

// Derive the block subtitle from the first three populated inputs, joined by
// " - " (mirroring Amplicon Alignment). Falls back to "Tite-Seq Analysis" when
// nothing is selected yet so the shelf label is never blank.
watchEffect(() => {
  const parts: string[] = [];

  const abundanceLabel = app.model.outputs.abundanceOptions?.find(
    (o) => app.model.data.abundanceRef && plRefsEqual(o.ref, app.model.data.abundanceRef),
  )?.label;
  if (abundanceLabel) parts.push(abundanceLabel);

  const concentrationLabel = app.model.outputs.concentrationOptions?.find(
    (o) =>
      app.model.data.concentrationColumnRef &&
      plRefsEqual(o.ref, app.model.data.concentrationColumnRef),
  )?.label;
  if (concentrationLabel) parts.push(concentrationLabel);

  const binLabel = app.model.outputs.binOptions?.find(
    (o) => app.model.data.binColumnRef && plRefsEqual(o.ref, app.model.data.binColumnRef),
  )?.label;
  if (binLabel) parts.push(binLabel);

  if (parts.length < 3 && app.model.data.targetAntigen) {
    parts.push(app.model.data.targetAntigen);
  }

  app.model.data.defaultBlockLabel = parts.slice(0, 3).join(" - ") || "Tite-Seq Analysis";
});
</script>

<template>
  <PlBtnGhost @click.stop="logOpen = true">
    Logs
    <template #append>
      <PlMaskIcon24 name="file-logs" />
    </template>
  </PlBtnGhost>
  <PlBtnGhost @click.stop="app.model.data.settingsOpen = true">
    Inputs
    <template #append>
      <PlMaskIcon24 name="settings" />
    </template>
  </PlBtnGhost>

  <SettingsDrawer v-model="app.model.data.settingsOpen" />
  <PlSlideModal v-model="logOpen" width="80%">
    <template #title>Fit log</template>
    <PlLogView :log-handle="app.model.outputs.logHandle" />
  </PlSlideModal>
</template>
