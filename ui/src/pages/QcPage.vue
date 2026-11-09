<script setup lang="ts">
/**
 * One page for every QC view, a tab each, as feature-integration's Panel QC page holds its views:
 * the section list stays at the two tables and this page.
 */
import type { QcView } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { qcViewsFor } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { PlTabs } from "@platforma-sdk/ui-vue";
import { computed, ref } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";
import AffinityVsFitView from "../qc/AffinityVsFitView.vue";
import KDDistributionView from "../qc/KDDistributionView.vue";
import ReadDepthView from "../qc/ReadDepthView.vue";
import ReplicateNoiseView from "../qc/ReplicateNoiseView.vue";
import SynonymousNoiseView from "../qc/SynonymousNoiseView.vue";
import TitrationCurvesView from "../qc/TitrationCurvesView.vue";
import TitrationEdgesView from "../qc/TitrationEdgesView.vue";
import VsParentView from "../qc/VsParentView.vue";

const app = useApp();

const VIEWS = {
  "vs-parent": VsParentView,
  "synonymous-noise": SynonymousNoiseView,
  "replicate-noise": ReplicateNoiseView,
  "titration-edges": TitrationEdgesView,
  "read-depth": ReadDepthView,
  "titration-curves": TitrationCurvesView,
  "kd-distribution": KDDistributionView,
  "affinity-vs-fit": AffinityVsFitView,
} satisfies Record<QcView, unknown>;

// While a run is in flight the run shape is unknown; the views every run has stay offered.
const tabs = computed(() => qcViewsFor(app.model.outputs.runShape ?? { synonymous: false }));

// Local, not stored. The chosen view is kept while a re-run hides its tab, and shown again when the
// tab returns; meanwhile the first tab shows.
const chosen = ref<QcView | undefined>(undefined);
const active = computed({
  get: (): QcView => tabs.value.find((t) => t.value === chosen.value)?.value ?? tabs.value[0].value,
  set: (view: QcView) => {
    chosen.value = view;
  },
});
</script>

<template>
  <TiteseqPage title="QC" mode="graph">
    <PlTabs v-model="active" :options="tabs" class="titeseq-qc-tabs" />
    <div class="titeseq-qc-view">
      <component :is="VIEWS[active]" />
    </div>
  </TiteseqPage>
</template>

<style scoped>
.titeseq-qc-tabs {
  margin: 8px 16px 0;
}
/* A plot fills its container, and a container with no height collapses to nothing: the view takes
   the height the tabs leave. */
.titeseq-qc-view {
  flex: 1;
  min-height: 480px;
  display: flex;
  flex-direction: column;
}
</style>
