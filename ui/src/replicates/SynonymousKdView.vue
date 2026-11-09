<script setup lang="ts">
/**
 * Synonymous mode with copies and replicates: one amino-acid variant per tab, each nt copy's
 * -log10 Kd,app as the mean over replicates with a standard-error bar. A replicate at an edge of
 * the titration contributes that edge.
 */
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import { SYNONYMOUS_KD_DEFAULT_STATE } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

const app = useApp();

const REPLICATE_AXIS = "pl7.app/vdj/replicate";

/**
 * The state must exist in data before the binding reads it: a getter that returns a new default
 * object per read never lets the chart settle, and it renders blank.
 */
app.model.data.graphStateSynonymousKd ??= { ...SYNONYMOUS_KD_DEFAULT_STATE };

const graphState = computed({
  get: (): GraphMakerState => app.model.data.graphStateSynonymousKd ?? SYNONYMOUS_KD_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data.graphStateSynonymousKd = state;
  },
});

const defaultOptions = computed((): PredefinedGraphOption<"discrete">[] | undefined => {
  const cols = app.model.outputs.synonymousKdPfCols;
  const kd = cols?.find((p) => p.spec.name === "pl7.app/vdj/kdNegLog10");
  // Each nt copy's amino-acid variant, on the nt axis: the chart cannot tab by the amino-acid
  // axis itself, which it reaches only through the linker, and falls back to the nt axis.
  const aaVariant = cols?.find((p) => p.spec.name === "pl7.app/vdj/aaVariantId");
  const ntAxis = kd?.spec.axesSpec.find((a) => a.name !== REPLICATE_AXIS);
  if (!kd || !aaVariant || !ntAxis) return undefined;
  return [
    { inputName: "y", selectedSource: kd.spec },
    { inputName: "primaryGrouping", selectedSource: ntAxis },
    { inputName: "tabBy", selectedSource: aaVariant.spec },
  ];
});
</script>

<template>
  <GraphMaker
    v-model="graphState"
    chart-type="discrete"
    :p-frame="app.model.outputs.synonymousKdPf"
    :default-options="defaultOptions"
    default-palette="bright"
    :status-text="{
      noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
    }"
  />
</template>
