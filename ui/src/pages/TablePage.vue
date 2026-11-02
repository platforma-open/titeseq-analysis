<script setup lang="ts">
import type { PTableKey } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { createPlDataTableStateV2 } from "@platforma-sdk/model";
import {
  PL_PLACEHOLDER_TEXTS,
  PlAgDataTableV2,
  PlPlaceholder,
  PlSlideModal,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.summaryTable,
});

// Replicate drill-down. The selected key lives in block data, outside args: the model filters the
// per-replicate frame to it. Written only on a user gesture.
const drilldownOpen = computed({
  get: () => app.model.data.replicateDrilldownKey !== undefined,
  set: (open: boolean) => {
    if (!open) app.model.data.replicateDrilldownKey = undefined;
  },
});

function openDrilldown(key?: PTableKey) {
  if (key === undefined || app.model.outputs.replicates !== true) return;
  app.model.data.replicateDrilldownTableState = undefined;
  app.model.data.replicateDrilldownKey = key;
}

const drilldownTableState = computed({
  get: () => app.model.data.replicateDrilldownTableState ?? createPlDataTableStateV2(),
  set: (value) => {
    app.model.data.replicateDrilldownTableState = value;
  },
});

// The source changes on a click, not on a run, so `sourceId` is required.
const drilldownSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.replicateDrilldownTable,
  sourceId: () => app.model.data.replicateDrilldownKey?.join(" "),
});

// Mounted only once the output holds a table (or an error): a grid mounted before that stays on
// "Data is not computed".
const drilldownReady = computed(() => {
  const out = app.model.outputs.replicateDrilldownTable;
  return out !== undefined && (!out.ok || out.value !== undefined);
});

const cellButtonAxis = computed(() =>
  app.model.outputs.replicates === true ? app.model.outputs.summaryKeyAxisId : undefined,
);
</script>

<template>
  <TiteseqPage title="Clonotype Fit Results" mode="table" main>
    <PlAgDataTableV2
      v-model="app.model.data.tableState"
      :settings="tableSettings"
      show-export-button
      :show-cell-button-for-axis-id="cellButtonAxis"
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      no-rows-text="No clonotypes available."
      @row-double-clicked="openDrilldown"
      @cell-button-clicked="openDrilldown"
    />
    <PlSlideModal v-model="drilldownOpen" width="80%">
      <!-- The variant is named by the label column inside the table below: the model cannot
           resolve it, since the axis value is an opaque key. -->
      <template #title>Replicates</template>
      <PlAgDataTableV2
        v-if="drilldownReady"
        v-model="drilldownTableState"
        :settings="drilldownSettings"
        show-export-button
        no-rows-text="No per-replicate values for this variant."
      />
      <PlPlaceholder
        v-else
        variant="table"
        :title="PL_PLACEHOLDER_TEXTS.LOADING.title"
        :subtitle="PL_PLACEHOLDER_TEXTS.LOADING.subtitle"
      />
    </PlSlideModal>
  </TiteseqPage>
</template>
