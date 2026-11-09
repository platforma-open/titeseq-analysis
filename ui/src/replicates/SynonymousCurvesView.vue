<script setup lang="ts">
/**
 * Synonymous mode with copies: one amino-acid variant per tab, each nt copy's fitted curve at every
 * concentration as a line through the mean over replicates, with a standard-error bar. A category
 * chart: the only chart that draws error bars, so the concentrations sit evenly spaced.
 */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import { SYNONYMOUS_CURVES_DEFAULT_STATE } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";
import { useConcentrationOrder } from "../composables/useConcentrationOrder";

const app = useApp();

const CONCENTRATION_AXIS = "pl7.app/vdj/concentration";
const REPLICATE_AXIS = "pl7.app/vdj/replicate";

const fitted = computed(() =>
  app.model.outputs.titrationCurvesPfCols?.find((p) => p.spec.name === "pl7.app/vdj/fittedMeanBin"),
);

const concentrationAxis = computed(() =>
  fitted.value?.spec.axesSpec.find((a) => a.name === CONCENTRATION_AXIS),
);

/**
 * The state must exist in data before the binding reads it: a getter that returns a new default
 * object per read never lets the chart settle, and it renders blank.
 */
app.model.data.graphStateSynonymousCurves ??= { ...SYNONYMOUS_CURVES_DEFAULT_STATE };

useConcentrationOrder(
  computed({
    get: () => app.model.data.graphStateSynonymousCurves,
    set: (state) => {
      app.model.data.graphStateSynonymousCurves = state;
    },
  }),
  concentrationAxis,
);

const graphState = computed({
  get: (): GraphMakerState =>
    app.model.data.graphStateSynonymousCurves ?? SYNONYMOUS_CURVES_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data.graphStateSynonymousCurves = state;
  },
});

const defaultOptions = computed((): PredefinedGraphOption<"discrete">[] | undefined => {
  const curveCols = app.model.outputs.titrationCurvesPfCols;
  // Each nt copy's amino-acid variant, on the nt axis: the chart cannot tab by the amino-acid
  // axis itself, which it reaches only through the linker, and falls back to the nt axis.
  const aaVariant = curveCols?.find((p) => p.spec.name === "pl7.app/vdj/aaVariantId");
  if (!fitted.value || !aaVariant) return undefined;
  const concentration = concentrationAxis.value;
  const ntAxis = fitted.value.spec.axesSpec.find(
    (a) => a.name !== CONCENTRATION_AXIS && a.name !== REPLICATE_AXIS,
  );
  if (!concentration || !ntAxis) return undefined;
  return [
    { inputName: "y", selectedSource: fitted.value.spec },
    { inputName: "primaryGrouping", selectedSource: concentration },
    { inputName: "secondaryGrouping", selectedSource: ntAxis },
    { inputName: "tabBy", selectedSource: aaVariant.spec },
  ];
});
</script>

<template>
  <GraphMaker
    v-model="graphState"
    chart-type="discrete"
    :data-state-key="app.model.outputs.titrationCurvesPf"
    :p-frame="app.model.outputs.titrationCurvesPf"
    :default-options="defaultOptions"
    default-palette="bright"
    :status-text="{
      noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
    }"
  />
</template>
