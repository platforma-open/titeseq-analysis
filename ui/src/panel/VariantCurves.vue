<script setup lang="ts">
/**
 * One variant's titration curves inside the variant panel: each replicate's observed points and
 * fitted curve on a log concentration axis. The open variant reaches the chart as a hidden pin
 * (`fixedOptions`), so one saved chart state serves every variant and opening one writes nothing.
 */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import {
  VARIANT_CURVES_DEFAULT_STATE,
  variantCurveDefaults,
  variantCurvePins,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

const props = defineProps<{ rowKey: string }>();

const app = useApp();

const cols = computed(() => app.model.outputs.titrationCurvesPfCols);

// The state must exist in data before the binding reads it, or the chart never settles.
app.model.data.graphStateVariantCurves ??= { ...VARIANT_CURVES_DEFAULT_STATE };

const graphState = computed({
  get: (): GraphMakerState =>
    app.model.data.graphStateVariantCurves ?? VARIANT_CURVES_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data.graphStateVariantCurves = state;
  },
});

const defaultOptions = computed((): PredefinedGraphOption<"scatterplot">[] | undefined =>
  variantCurveDefaults(cols.value?.map((c) => c.spec) ?? []),
);

const pins = computed(() =>
  variantCurvePins(
    cols.value?.map((c) => c.spec) ?? [],
    props.rowKey,
    app.model.outputs.synonymous === true,
  ),
);
</script>

<template>
  <div class="titeseq-variant-curves">
    <GraphMaker
      v-if="pins"
      :key="rowKey"
      v-model="graphState"
      chart-type="scatterplot"
      :p-frame="app.model.outputs.titrationCurvesPf"
      :default-options="defaultOptions"
      :fixed-options="pins"
      default-palette="bright"
    />
    <div v-else class="titeseq-variant-curves-none">No titration curves for this variant.</div>
  </div>
</template>

<style scoped>
.titeseq-variant-curves {
  flex: 1;
  display: flex;
  flex-direction: column;
  min-height: 480px;
}
.titeseq-variant-curves-none {
  padding: 16px;
  color: var(--txt-03);
}
</style>
