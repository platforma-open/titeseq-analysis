<script setup lang="ts">
import { PlBtnGhost, PlLogView, PlMaskIcon24, PlSlideModal } from "@platforma-sdk/ui-vue";
import { onMounted, ref } from "vue";
import { useApp } from "../app";
import SettingsDrawer from "./SettingsDrawer.vue";

const app = useApp();
const logOpen = ref(false);

onMounted(() => {
  if (app.model.data.abundanceRef === undefined && !app.model.data.settingsOpen) {
    app.model.data.settingsOpen = true;
  }
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
