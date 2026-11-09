<script setup lang="ts">
/**
 * Runs with sample replicates: every value each replicate contributed, a tab per view, as the QC
 * page holds its views.
 */
import type { ReplicateView } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { replicateViewsFor } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { PlTabs } from "@platforma-sdk/ui-vue";
import { computed, ref } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";
import NtReplicatesView from "../replicates/NtReplicatesView.vue";
import SynonymousCurvesView from "../replicates/SynonymousCurvesView.vue";
import SynonymousKdView from "../replicates/SynonymousKdView.vue";
import VariantReplicatesView from "../replicates/VariantReplicatesView.vue";

const app = useApp();

const VIEWS = {
  "variant-replicates": VariantReplicatesView,
  "nt-replicates": NtReplicatesView,
  "synonymous-curves": SynonymousCurvesView,
  "synonymous-kd": SynonymousKdView,
} satisfies Record<ReplicateView, unknown>;

// While a run is in flight the run shape is unknown; the view every run has stays offered.
const tabs = computed(() => replicateViewsFor(app.model.outputs.runShape ?? { synonymous: false }));

// Local, not stored. The chosen view is kept while a re-run hides its tab.
const chosen = ref<ReplicateView | undefined>(undefined);
const active = computed({
  get: (): ReplicateView =>
    tabs.value.find((t) => t.value === chosen.value)?.value ?? tabs.value[0].value,
  set: (view: ReplicateView) => {
    chosen.value = view;
  },
});
</script>

<template>
  <TiteseqPage title="Replicates" mode="graph">
    <PlTabs v-model="active" :options="tabs" class="titeseq-replicate-tabs" />
    <div class="titeseq-replicate-view">
      <component :is="VIEWS[active]" />
    </div>
  </TiteseqPage>
</template>

<style scoped>
.titeseq-replicate-tabs {
  margin: 8px 16px 0;
}
/* A plot or grid fills its container, and a container with no height collapses to nothing: the
   view takes the height the tabs leave. */
.titeseq-replicate-view {
  flex: 1;
  min-height: 480px;
  display: flex;
  flex-direction: column;
}
</style>
