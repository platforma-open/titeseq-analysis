<script setup lang="ts">
import type { PTableKey } from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  createPlDataTableStateV2,
  getRawPlatformaInstance,
  isLabelColumn,
  pTableValue,
} from "@platforma-sdk/model";
import type { PTableHandle } from "@platforma-sdk/model";
import {
  PL_PLACEHOLDER_TEXTS,
  PlAgDataTableV2,
  PlPlaceholder,
  PlSlideModal,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
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

// The panel title names the variant. The model cannot resolve the label, so the UI reads it from
// the drill-down table's hidden label column: the table holds one variant, so row 0 names it.
async function readVariantLabel(handle: PTableHandle): Promise<string | undefined> {
  const driver = getRawPlatformaInstance().pFrameDriver;
  const specs = await driver.getSpec(handle);
  const iLabel = specs.findIndex((s) => s.type === "column" && isLabelColumn(s.spec));
  if (iLabel < 0) return undefined;
  const { rows } = await driver.getShape(handle);
  if (rows === 0) return undefined;
  const [labels] = await driver.getData(handle, [iLabel], { offset: 0, length: 1 });
  const label = pTableValue(labels, 0);
  return typeof label === "string" ? label : undefined;
}

const drilldownVariantLabel = ref<string>();
watch(
  () => {
    const out = app.model.outputs.replicateDrilldownTable;
    return out?.ok ? out.value?.fullTableHandle : undefined;
  },
  async (handle) => {
    drilldownVariantLabel.value = undefined;
    if (handle === undefined) return;
    const label = await readVariantLabel(handle);
    // A later selection may have replaced the table while the read was in flight.
    const current = app.model.outputs.replicateDrilldownTable;
    if (current?.ok && current.value?.fullTableHandle === handle) {
      drilldownVariantLabel.value = label;
    }
  },
  { immediate: true },
);

const cellButtonAxis = computed(() =>
  app.model.outputs.replicates === true ? app.model.outputs.summaryKeyAxisId : undefined,
);
// PlAgDataTableV2 reads `show-cell-button-for-axis-id` only when its table data changes. The two
// outputs behind the axis can resolve after the table, so the key remounts the grid when the axis
// arrives. Without it the Open button can be missing.
const cellButtonKey = computed(() => JSON.stringify(cellButtonAxis.value ?? null));
</script>

<template>
  <TiteseqPage title="Clonotype Fit Results" mode="table" main>
    <PlAgDataTableV2
      :key="cellButtonKey"
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
      <template #title>{{ drilldownVariantLabel ?? "Replicates" }}</template>
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
