<script setup lang="ts">
import type { PTableKey } from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  createPlDataTableStateV2,
  getRawPlatformaInstance,
  isLabelColumn,
  pTableValue,
} from "@platforma-sdk/model";
import type { PTableColumnSpec, PTableHandle } from "@platforma-sdk/model";
import {
  PL_PLACEHOLDER_TEXTS,
  PlAgDataTableV2,
  PlAgTextAndButtonCell,
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

// The Open button on the variant column. A custom renderer rather than `show-cell-button-for-axis-id`:
// PlAgDataTableV2 reads that prop only when its table data changes, so an axis id that resolves after
// the table never gets a button. `cell-renderer-selector` is re-applied whenever it changes.
// `invokeRowsOnDoubleClick` makes the button fire the row's double-click, so the button and a
// double-click anywhere on the row take the same path into the drill-down.
const keyAxisName = computed(() =>
  app.model.outputs.replicates === true ? app.model.outputs.summaryKeyAxisId?.name : undefined,
);

function isVariantColumn(spec: PTableColumnSpec | undefined, axisName: string): boolean {
  if (spec === undefined) return false;
  if (spec.type === "axis") return spec.spec.name === axisName;
  return (
    isLabelColumn(spec.spec) &&
    spec.spec.axesSpec.length === 1 &&
    spec.spec.axesSpec[0].name === axisName
  );
}

const cellRendererSelector = computed(() => {
  const axisName = keyAxisName.value;
  if (axisName === undefined) return undefined;
  return (params: { colDef?: { context?: PTableColumnSpec } }) =>
    isVariantColumn(params.colDef?.context, axisName)
      ? { component: PlAgTextAndButtonCell, params: { invokeRowsOnDoubleClick: true } }
      : undefined;
});
</script>

<template>
  <TiteseqPage title="Clonotype Fit Results" mode="table" main>
    <p v-if="keyAxisName !== undefined" class="drilldown-hint">
      Double-click a variant, or click Open, to see its replicates.
    </p>
    <PlAgDataTableV2
      v-model="app.model.data.tableState"
      :settings="tableSettings"
      show-export-button
      :cell-renderer-selector="cellRendererSelector"
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      no-rows-text="No clonotypes available."
      @row-double-clicked="openDrilldown"
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

<style scoped>
.drilldown-hint {
  margin: 0 0 8px;
  color: var(--txt-03);
  font-size: 13px;
}
</style>
