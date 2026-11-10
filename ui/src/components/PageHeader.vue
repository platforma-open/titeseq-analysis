<script setup lang="ts">
import type { ImportFileHandle } from "@platforma-sdk/model";
import type { FileExportEntry } from "@platforma-sdk/ui-vue";
import {
  PlBtnExportArchive,
  PlBtnGhost,
  PlLogView,
  PlMaskIcon24,
  PlSlideModal,
} from "@platforma-sdk/ui-vue";
import { computed, onMounted, ref } from "vue";
import { useApp } from "../app";
import SettingsDrawer from "./SettingsDrawer.vue";

const props = defineProps<{ showSettings?: boolean }>();

const app = useApp();
const logOpen = ref(false);

// Every table the run wrote, nucleotide variants and replicates included, as the software wrote it.
const rawFileExports = computed<FileExportEntry[]>(() =>
  (app.model.outputs.rawFiles ?? []).flatMap(({ fileName, handle }) =>
    handle === undefined
      ? []
      : [{ importHandle: fileName as ImportFileHandle, blobHandle: handle, fileName }],
  ),
);

onMounted(() => {
  if (
    props.showSettings &&
    app.model.data.abundanceRef === undefined &&
    !app.model.data.settingsOpen
  ) {
    app.model.data.settingsOpen = true;
  }
});
</script>

<template>
  <template v-if="props.showSettings">
    <PlBtnExportArchive
      :file-exports="rawFileExports"
      suggested-file-name="titeseq-raw-data.zip"
      :disabled="rawFileExports.length === 0"
      :file-picker-types="[{ description: 'ZIP files', accept: { 'application/zip': ['.zip'] } }]"
    >
      Export raw data
    </PlBtnExportArchive>
    <PlBtnGhost @click.stop="logOpen = true">
      Logs
      <template #append>
        <PlMaskIcon24 name="file-logs" />
      </template>
    </PlBtnGhost>
    <PlBtnGhost @click.stop="app.model.data.settingsOpen = true">
      Settings
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
