<script setup lang="ts">
import { PL_PLACEHOLDER_TEXTS, PlAlert, PlPlaceholder } from "@platforma-sdk/ui-vue";
import { computed, ref } from "vue";
import { useApp } from "../app";

// The Run Quality page reads counts the run writes. A run made before the block wrote them has
// results but no counts, and says so.
const app = useApp();

const state = computed<"ready" | "running" | "noRun" | "oldRun">(() => {
  if (app.model.outputs.runQuality !== undefined) return "ready";
  if (app.model.outputs.isRunning) return "running";
  return app.model.outputs.isEmpty === undefined ? "noRun" : "oldRun";
});

// Without a known-variant list the counts cannot tell a designed variant from a sequencing error.
const leftOut = computed(() => {
  const quality = app.model.outputs.runQuality;
  return quality?.scope === "attempted" ? quality.leftOutReasons.join(", ") : undefined;
});
// Closed until the block view is opened again, like the settings warnings.
const leftOutShown = ref(true);
</script>

<template>
  <template v-if="state === 'ready'">
    <PlAlert
      v-if="leftOut !== undefined"
      v-model="leftOutShown"
      type="warn"
      closeable
      class="run-quality-gate run-quality-gate--closeable"
    >
      The input dataset has no known-variant list, so these graphs cannot tell the variants you
      designed from sequencing errors. They count the variants with enough data for a fit to be
      attempted, and leave out the fits that failed for: {{ leftOut }}. The Raw Counts tab lists
      every variant. To count only your known variants, add a known-variant list in the block that
      produced the dataset.
    </PlAlert>
    <slot />
  </template>
  <PlPlaceholder
    v-else-if="state === 'running'"
    variant="graph"
    :title="PL_PLACEHOLDER_TEXTS.RUNNING.title"
    :subtitle="PL_PLACEHOLDER_TEXTS.RUNNING.subtitle"
  />
  <PlAlert v-else-if="state === 'noRun'" type="info" class="run-quality-gate">
    Configure the block in Settings on the Main page and run it.
  </PlAlert>
  <PlAlert v-else type="warn" class="run-quality-gate">
    These results were made before the block counted fit outcomes. Run the block again to see them.
  </PlAlert>
</template>

<style scoped>
.run-quality-gate {
  margin: 16px;
}
/* Keeps the text clear of the close button. */
.run-quality-gate--closeable {
  padding-right: 48px;
}
</style>
