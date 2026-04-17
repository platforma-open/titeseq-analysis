<script setup lang="ts">
import { PlBtnGhost, PlLogView, PlMaskIcon24, PlSlideModal } from "@platforma-sdk/ui-vue";
import { ref } from "vue";
import { useApp } from "../app";
import SettingsDrawer from "./SettingsDrawer.vue";

const app = useApp();
const settingsOpen = ref(app.model.args.abundanceRef === undefined);
const logOpen = ref(false);
</script>

<template>
  <PlBtnGhost @click.stop="logOpen = true">
    Logs
    <template #append>
      <PlMaskIcon24 name="file-logs" />
    </template>
  </PlBtnGhost>
  <PlBtnGhost @click.stop="settingsOpen = true">
    Settings
    <template #append>
      <PlMaskIcon24 name="settings" />
    </template>
  </PlBtnGhost>

  <SettingsDrawer v-model="settingsOpen" />
  <PlSlideModal v-model="logOpen" width="80%">
    <template #title>Fit log</template>
    <PlLogView :log-handle="app.model.outputs.logHandle" />
  </PlSlideModal>
</template>
