<script setup lang="ts">
import type { PTableKey } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { drilldownShownColumns } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { createPlDataTableStateV2, isLabelColumn } from "@platforma-sdk/model";
import type { PTableColumnSpec } from "@platforma-sdk/model";
import {
  PL_PLACEHOLDER_TEXTS,
  PlAgDataTableV2,
  PlAgTextAndButtonCell,
  PlDropdownMulti,
  PlPlaceholder,
  PlSlideModal,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";
import { useNoun } from "../composables/useNoun";
import { readDrilldownFacts } from "../utils/readDrilldownFacts";

const app = useApp();
const { noun, Noun } = useNoun();

// Synonymous mode shows one row per amino-acid variant.
const synonymous = computed(() => app.model.outputs.runShape?.synonymous === true);
const replicates = computed(() => app.model.outputs.runShape?.replicates === true);
const title = computed(() =>
  synonymous.value ? "Amino Acid Variant Fit Results" : `${Noun.value} Fit Results`,
);
const noRowsText = computed(() =>
  synonymous.value ? "No amino acid variants available." : `No ${noun.value}s available.`,
);

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
  if (key === undefined || !replicates.value) return;
  app.model.data.replicateDrilldownTableState = undefined;
  app.model.data.replicateDrilldownEmptyColumns = undefined;
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
// "Data is not computed". It also waits until the table hides the columns the panel's list leaves
// out: a grid mounted on an earlier table stores that table's columns as its own state, and
// keeps them.
const drilldownReady = computed(() => {
  const out = app.model.outputs.replicateDrilldownTable;
  if (out === undefined) return false;
  if (!out.ok) return true;
  if (out.value === undefined || app.model.data.replicateDrilldownEmptyColumns === undefined)
    return false;
  const hiddenInTable = Object.values(out.value.columnsMeta?.columns ?? {}).filter(
    (c) => c.visibility === "optional",
  ).length;
  return hiddenInTable === drilldownColumnOptions.value.length - drilldownShown.value.length;
});

// The panel title names the variant, read from the drill-down table's hidden label column. The
// same read finds the columns with no value for this variant, which the model then leaves out of
// the default columns.
const drilldownVariantLabel = ref<string>();
watch(
  () => {
    const out = app.model.outputs.replicateDrilldownTable;
    return out?.ok ? out.value?.fullTableHandle : undefined;
  },
  async (handle) => {
    drilldownVariantLabel.value = undefined;
    if (handle === undefined) return;
    const { label, emptyColumns } = await readDrilldownFacts(handle);
    // A later selection may have replaced the table while the read was in flight.
    const current = app.model.outputs.replicateDrilldownTable;
    if (current?.ok && current.value?.fullTableHandle === handle) {
      drilldownVariantLabel.value = label;
      const known = app.model.data.replicateDrilldownEmptyColumns;
      if (known === undefined || known.join("\n") !== emptyColumns.join("\n")) {
        app.model.data.replicateDrilldownEmptyColumns = emptyColumns;
      }
    }
  },
  { immediate: true },
);

// The panel's own column list: the grid's column button sits in the page title, behind the panel.
const drilldownColumnOptions = computed(() =>
  (app.model.outputs.replicateDrilldownColumns ?? []).map((c) => ({
    value: c.key,
    label: c.label,
  })),
);
const drilldownShown = computed({
  get: () =>
    drilldownShownColumns(
      drilldownColumnOptions.value.map((o) => o.value),
      app.model.data.replicateDrilldownEmptyColumns,
      app.model.data.replicateDrilldownColumnChoice,
    ),
  set: (shown: string[]) => {
    const before = drilldownShown.value;
    const choice = { ...app.model.data.replicateDrilldownColumnChoice };
    for (const { value } of drilldownColumnOptions.value) {
      if (shown.includes(value) !== before.includes(value)) choice[value] = shown.includes(value);
    }
    app.model.data.replicateDrilldownColumnChoice = choice;
    // The grid keeps its own hidden columns once it has a state; without one it takes the model's.
    app.model.data.replicateDrilldownTableState = undefined;
  },
});
// A new set of default columns remounts the grid, which reads the defaults on mount only.
const drilldownGridKey = computed(() => drilldownShown.value.join("\n"));

// The Open button on the variant column. A custom renderer rather than `show-cell-button-for-axis-id`:
// PlAgDataTableV2 reads that prop only when its table data changes, so an axis id that resolves after
// the table never gets a button. `cell-renderer-selector` is re-applied whenever it changes.
// `invokeRowsOnDoubleClick` makes the button fire the row's double-click, so the button and a
// double-click anywhere on the row take the same path into the drill-down.
const keyAxisName = computed(() =>
  replicates.value ? app.model.outputs.summaryKeyAxisId?.name : undefined,
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
  <TiteseqPage :title="title" mode="table" main>
    <PlAgDataTableV2
      v-model="app.model.data.tableState"
      :settings="tableSettings"
      show-export-button
      :cell-renderer-selector="cellRendererSelector"
      not-ready-text="Open Settings (top right) to configure the block and run it."
      :no-rows-text="noRowsText"
      @row-double-clicked="openDrilldown"
    />
    <PlSlideModal v-model="drilldownOpen" width="80%">
      <template #title>{{ drilldownVariantLabel ?? "Replicates" }}</template>
      <PlDropdownMulti
        v-if="drilldownColumnOptions.length > 0"
        v-model="drilldownShown"
        :options="drilldownColumnOptions"
        label="Columns"
      />
      <PlAgDataTableV2
        v-if="drilldownReady"
        :key="drilldownGridKey"
        v-model="drilldownTableState"
        :settings="drilldownSettings"
        disable-columns-panel
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
