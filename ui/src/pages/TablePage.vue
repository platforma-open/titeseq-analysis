<script setup lang="ts">
import type { PTableKey } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { canonicalizeJson, createPlDataTableStateV2, isLabelColumn } from "@platforma-sdk/model";
import type { PTableColumnSpec } from "@platforma-sdk/model";
import {
  PL_PLACEHOLDER_TEXTS,
  PlAgDataTableV2,
  PlAgTextAndButtonCell,
  PlPlaceholder,
  PlSlideModal,
  PlTabs,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";
import { useNoun } from "../composables/useNoun";
import VariantCurves from "../panel/VariantCurves.vue";
import { handleOf, loadCurveVariantLabel } from "../panel/panelVariants";
import { readDefaultHiddenColumns, seededTableState } from "../utils/emptyColumns";
import { readTableLabel } from "../utils/readTableLabel";

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
  if (key === undefined) return;
  app.model.data.replicateDrilldownTableState = undefined;
  app.model.data.replicateDrilldownKey = key;
  // Only the client that opened the row seeds its table. The settings canonicalize the source id.
  seedFor = canonicalizeJson(key.join(" "));
  seededSource.value = undefined;
}

// The variant panel's tabs. Curves exist for every run, Replicates only with replicates. A single
// tab shows none.
type PanelTab = "curves" | "replicates";
const panelTabs = computed(
  () =>
    [
      { value: "curves" as const, label: "Curves" },
      ...(replicates.value ? [{ value: "replicates" as const, label: "Replicates" }] : []),
    ] satisfies { value: PanelTab; label: string }[],
);
const chosenPanelTab = ref<PanelTab>("curves");
// A row without curves opens on its replicates, until a tab is picked for it.
const curvesEmpty = ref(false);
const pickedForRow = ref(false);
const panelTab = computed({
  get: (): PanelTab => {
    const has = (tab: PanelTab) => panelTabs.value.some((t) => t.value === tab);
    if (curvesEmpty.value && !pickedForRow.value && has("replicates")) return "replicates";
    return has(chosenPanelTab.value) ? chosenPanelTab.value : "curves";
  },
  set: (tab: PanelTab) => {
    chosenPanelTab.value = tab;
    pickedForRow.value = true;
  },
});
const openRowKey = computed(() => {
  const value = app.model.data.replicateDrilldownKey?.[0];
  return value === undefined || value === null ? undefined : String(value);
});
watch(openRowKey, () => {
  curvesEmpty.value = false;
  pickedForRow.value = false;
});

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

// The source whose grid may mount: its columns empty in every replicate are seeded as hidden, or
// it needs no seeding.
const seededSource = ref<string>();
// The source this client opened and has still to seed. A client that did not open the row seeds
// nothing: it shows the state the opening client wrote, so two clients never seed one variant.
let seedFor: string | undefined;
// Watches the handle, not the output: the output object is updated in place when the next
// variant's table arrives, so a watch on it does not fire again.
//
// An output -> data write that cannot loop. It writes `replicateDrilldownTableState`, which
// `replicateDrilldownTable` reads, under three guards:
// - A user gesture starts it. Only `openDrilldown` sets `seedFor`, a local variable, so a client
//   that did not open the row never writes.
// - It writes once per gesture. `seedFor` is cleared in the same tick as the write, so the table
//   the write changes cannot start a second write.
// - The state it writes is UI state. It reaches this table only, never args or the trace, so it
//   cannot stale the block.
watch(
  [
    () => drilldownSettings.value.sourceId,
    () => {
      const out = app.model.outputs.replicateDrilldownTable;
      if (out === undefined) return undefined;
      return out.ok ? out.value?.fullTableHandle : "error";
    },
  ],
  async ([sourceId, handle], _old, onCleanup) => {
    let stale = false;
    onCleanup(() => {
      stale = true;
    });
    if (sourceId === null || handle === undefined) return;
    if (seedFor !== sourceId) {
      seededSource.value = sourceId;
      return;
    }
    const out = app.model.outputs.replicateDrilldownTable;
    if (handle === "error" || !out?.ok) {
      // The output can still hold the previous variant's error. The grid shows it; the row is
      // seeded when its table arrives.
      seededSource.value = sourceId;
      return;
    }
    const key = openRowKey.value;
    if (key === undefined) return;
    // Unmounts a grid shown on the error, so it remounts on the seeded state.
    seededSource.value = undefined;
    let hidden;
    try {
      hidden = await readDefaultHiddenColumns(handle, key, out.value?.columnsMeta);
    } catch (e) {
      // The table still opens, without the empty columns hidden.
      if (!stale) {
        seedFor = undefined;
        seededSource.value = sourceId;
      }
      throw e;
    }
    if (stale || hidden === undefined) return;
    if (hidden.length > 0) {
      app.model.data.replicateDrilldownTableState = seededTableState(sourceId, hidden);
    }
    seedFor = undefined;
    seededSource.value = sourceId;
  },
  { immediate: true },
);

// Mounted only once the output holds a table (or an error): a grid mounted before that stays on
// "Data is not computed". It also waits for the seeded state: a grid mounted first saves its own.
const drilldownReady = computed(() => {
  const out = app.model.outputs.replicateDrilldownTable;
  return (
    out !== undefined &&
    (!out.ok || out.value !== undefined) &&
    seededSource.value === drilldownSettings.value.sourceId
  );
});

// The panel title names the variant, read from the drill-down table's hidden label column.
const drilldownVariantLabel = ref<string>();
watch(
  () => {
    const out = app.model.outputs.replicateDrilldownTable;
    return out?.ok ? out.value?.fullTableHandle : undefined;
  },
  async (handle) => {
    drilldownVariantLabel.value = undefined;
    if (handle === undefined) return;
    const label = await readTableLabel(handle);
    // A later selection may have replaced the table while the read was in flight.
    const current = app.model.outputs.replicateDrilldownTable;
    if (current?.ok && current.value?.fullTableHandle === handle) {
      drilldownVariantLabel.value = label;
    }
  },
  { immediate: true },
);

// Without replicates there is no drill-down table to read the label from. Outside synonymous mode
// the row is the curves' own variant, and the graph frame carries its label.
const curveVariantLabel = ref<string>();
watch(
  [openRowKey, () => handleOf(app.model.outputs.titrationCurvesPf), replicates, synonymous],
  async ([key, frame, replicates, synonymous], _old, onCleanup) => {
    let stale = false;
    onCleanup(() => {
      stale = true;
    });
    curveVariantLabel.value = undefined;
    if (key === undefined || frame === undefined || replicates || synonymous) return;
    const label = await loadCurveVariantLabel(frame, key);
    if (!stale) curveVariantLabel.value = label;
  },
  { immediate: true },
);
const panelLabel = computed(
  () => drilldownVariantLabel.value ?? curveVariantLabel.value ?? "Variant",
);

// The Open button on the variant column. A custom renderer rather than `show-cell-button-for-axis-id`:
// PlAgDataTableV2 reads that prop only when its table data changes, so an axis id that resolves after
// the table never gets a button. `cell-renderer-selector` is re-applied whenever it changes.
// `invokeRowsOnDoubleClick` makes the button fire the row's double-click, so the button and a
// double-click anywhere on the row take the same path into the drill-down.
const keyAxisName = computed(() => app.model.outputs.summaryKeyAxisId?.name);

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
      :cell-renderer-selector="cellRendererSelector"
      not-ready-text="Open Settings (top right) to configure the block and run it."
      :no-rows-text="noRowsText"
      @row-double-clicked="openDrilldown"
    />
    <PlSlideModal v-model="drilldownOpen" width="80%">
      <template #title>{{ panelLabel }}</template>
      <PlTabs v-if="panelTabs.length > 1" v-model="panelTab" :options="panelTabs" />
      <VariantCurves
        v-if="panelTab === 'curves' && openRowKey !== undefined"
        :row-key="openRowKey"
        @empty="curvesEmpty = $event"
      />
      <PlAgDataTableV2
        v-else-if="drilldownReady"
        v-model="drilldownTableState"
        :settings="drilldownSettings"
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
