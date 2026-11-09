<script setup lang="ts">
/**
 * Synonymous mode only. Three tabs. Variants: the amino-acid variants; double-click a row for its
 * nucleotide copies in a side panel, and with replicates, double-click a copy for its value in each
 * replicate. Curves (mean ± SE): each variant's observed and fitted values combined over
 * replicates, on evenly spaced concentrations. Curves (log axis): the observed mean and the mean
 * curve on a log concentration axis, without error bars.
 *
 * The open keys live in block data, outside args: the model filters each table to them. Written
 * only on a user gesture.
 */
import type { PTableKey } from "@platforma-open/platforma-open.titeseq-analysis.model";
import type { PlDataTableModel } from "@platforma-sdk/model";
import { createPlDataTableStateV2 } from "@platforma-sdk/model";
import {
  PL_PLACEHOLDER_TEXTS,
  PlAgDataTableV2,
  PlBtnSecondary,
  PlPlaceholder,
  PlSlideModal,
  PlTabs,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";
import CombinedCurvesView from "../explorer/CombinedCurvesView.vue";
import ScatterCurvesView from "../explorer/ScatterCurvesView.vue";
import { readTableLabel } from "../utils/readTableLabel";

const app = useApp();

const TABS = [
  { value: "variants" as const, label: "Variants" },
  { value: "curves" as const, label: "Curves (mean ± SE)" },
  { value: "scatter" as const, label: "Curves (log axis)" },
];
// Local, not stored: the open tab is view state.
const tab = ref<(typeof TABS)[number]["value"]>("variants");

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.explorerTable,
});

const panelOpen = computed({
  get: () => app.model.data.explorerAaKey !== undefined,
  set: (open: boolean) => {
    if (open) return;
    app.model.data.explorerAaKey = undefined;
    app.model.data.explorerNtKey = undefined;
  },
});

function openAa(key?: PTableKey) {
  if (key === undefined) return;
  app.model.data.explorerNtTableState = undefined;
  app.model.data.explorerNtKey = undefined;
  app.model.data.explorerAaKey = key;
}

function openNt(key?: PTableKey) {
  if (key === undefined || app.model.outputs.replicates !== true) return;
  app.model.data.explorerReplicateTableState = undefined;
  app.model.data.explorerNtKey = key;
}

function back() {
  app.model.data.explorerNtKey = undefined;
}

const ntOpen = computed(() => app.model.data.explorerNtKey !== undefined);

// The sources change on a click, not on a run, so `sourceId` is required.
const ntSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.explorerNtTable,
  sourceId: () => app.model.data.explorerAaKey?.join(" "),
});
const replicateSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.explorerReplicateTable,
  sourceId: () => app.model.data.explorerNtKey?.join(" "),
});

const tableState = computed({
  get: () => app.model.data.explorerTableState ?? createPlDataTableStateV2(),
  set: (value) => {
    app.model.data.explorerTableState = value;
  },
});
const ntTableState = computed({
  get: () => app.model.data.explorerNtTableState ?? createPlDataTableStateV2(),
  set: (value) => {
    app.model.data.explorerNtTableState = value;
  },
});
const replicateTableState = computed({
  get: () => app.model.data.explorerReplicateTableState ?? createPlDataTableStateV2(),
  set: (value) => {
    app.model.data.explorerReplicateTableState = value;
  },
});

type TableOutput = { ok: true; value?: PlDataTableModel } | { ok: false } | undefined;

// A grid mounted before its output holds a table stays on "Data is not computed".
const isReady = (out: TableOutput) => out !== undefined && (!out.ok || out.value !== undefined);
const ntReady = computed(() => isReady(app.model.outputs.explorerNtTable));
const replicateReady = computed(() => isReady(app.model.outputs.explorerReplicateTable));

/** The label of the variant a drill-down table holds, read from its hidden label column. */
function useTableLabel(output: () => TableOutput) {
  const label = ref<string>();
  watch(
    () => {
      const out = output();
      return out?.ok ? out.value?.fullTableHandle : undefined;
    },
    async (handle) => {
      label.value = undefined;
      if (handle === undefined) return;
      const read = await readTableLabel(handle);
      // A later selection may have replaced the table while the read was in flight.
      const current = output();
      if (current?.ok && current.value?.fullTableHandle === handle) label.value = read;
    },
    { immediate: true },
  );
  return label;
}

const ntLabel = useTableLabel(() => app.model.outputs.explorerReplicateTable);

const panelTitle = computed(() =>
  ntOpen.value ? `Replicates of ${ntLabel.value ?? "the copy"}` : "Nucleotide copies",
);
</script>

<template>
  <TiteseqPage title="Explorer" :mode="tab === 'variants' ? 'table' : 'graph'">
    <PlTabs v-model="tab" :options="TABS" class="titeseq-explorer-tabs" />
    <CombinedCurvesView v-if="tab === 'curves'" />
    <ScatterCurvesView v-else-if="tab === 'scatter'" />
    <PlAgDataTableV2
      v-else
      v-model="tableState"
      :settings="tableSettings"
      show-export-button
      not-ready-text="Open Inputs (top right) to configure the block and run it."
      no-rows-text="No amino acid variants available."
      @row-double-clicked="openAa"
    />
    <PlSlideModal v-model="panelOpen" width="80%">
      <template #title>{{ panelTitle }}</template>
      <template v-if="ntOpen">
        <PlBtnSecondary class="titeseq-explorer-back" @click="back">
          Back to nucleotide copies
        </PlBtnSecondary>
        <PlAgDataTableV2
          v-if="replicateReady"
          v-model="replicateTableState"
          :settings="replicateSettings"
          show-export-button
          no-rows-text="No per-replicate values for this copy."
        />
        <PlPlaceholder
          v-else
          variant="table"
          :title="PL_PLACEHOLDER_TEXTS.LOADING.title"
          :subtitle="PL_PLACEHOLDER_TEXTS.LOADING.subtitle"
        />
      </template>
      <template v-else>
        <PlAgDataTableV2
          v-if="ntReady"
          v-model="ntTableState"
          :settings="ntSettings"
          show-export-button
          no-rows-text="No nucleotide copies for this variant."
          @row-double-clicked="openNt"
        />
        <PlPlaceholder
          v-else
          variant="table"
          :title="PL_PLACEHOLDER_TEXTS.LOADING.title"
          :subtitle="PL_PLACEHOLDER_TEXTS.LOADING.subtitle"
        />
      </template>
    </PlSlideModal>
  </TiteseqPage>
</template>

<style scoped>
.titeseq-explorer-tabs {
  margin: 8px 16px 0;
}
.titeseq-explorer-back {
  align-self: flex-start;
  margin-bottom: 8px;
}
</style>
