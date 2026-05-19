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

// Auto-derive the subtitle placeholder. Differentiating field first so
// deriveLabels() downstream has a stable disambiguator:
//   1. targetAntigen — per-instance selector in multi-antigen studies
//   2. antigenColumn label — fallback when targetAntigen is unset
//   3. abundance label — distinguishes datasets / chains
// Bin and concentration columns are dropped: identical across instances on
// the same dataset, so they crowd out the differentiator without adding signal.
watchEffect(() => {
  const parts: string[] = [];

  if (app.model.data.targetAntigen) {
    parts.push(app.model.data.targetAntigen);
  } else if (app.model.data.antigenColumnRef) {
    const antigenLabel = app.model.outputs.antigenOptions?.find(
      (o) => app.model.data.antigenColumnRef && plRefsEqual(o.ref, app.model.data.antigenColumnRef),
    )?.label;
    if (antigenLabel) parts.push(antigenLabel);
  }

  const abundanceLabel = app.model.outputs.abundanceOptions?.find(
    (o) => app.model.data.abundanceRef && plRefsEqual(o.ref, app.model.data.abundanceRef),
  )?.label;
  if (abundanceLabel) parts.push(abundanceLabel);

  app.model.data.defaultBlockLabel = parts.join(" - ") || "Tite-Seq Analysis";
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

<style scoped>
.titeseq-facs-badge {
  display: inline-flex;
  align-items: center;
  padding: 2px 8px;
  border-radius: 10px;
  font-size: 11px;
  font-weight: 600;
  letter-spacing: 0.02em;
  background: var(--color-accent-subtle, rgba(46, 160, 67, 0.15));
  color: var(--color-accent, #2ea043);
  margin-right: 8px;
}
</style>
