<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const defaultOptions = computed((): PredefinedGraphOption<"scatterplot">[] | undefined => {
  const pCols = app.model.outputs.summaryPfCols;
  if (!pCols) return undefined;

  const kd = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/kd");
  const hillCoef = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/hillCoefficient",
  );
  const affinityClass = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/affinityClass",
  );
  const kdSe = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/kdLog10Se");
  const y = kdSe ?? hillCoef;
  if (!kd || !y || !affinityClass) return undefined;

  // y is the measure the equation that produced these results rates fit quality on: SE of
  // log10 Kd,app (present only under log_fluorescence), or the Hill coefficient. Read from the
  // results, not from `data.equation`, which runs ahead of them until the block is re-run. A fit
  // with no Kd,app, or no value on y, is not drawn; the table lists it with its failure reason.
  const options: PredefinedGraphOption<"scatterplot">[] = [
    { inputName: "x", selectedSource: kd.spec },
    { inputName: "y", selectedSource: y.spec },
    { inputName: "grouping", selectedSource: affinityClass.spec },
  ];
  // GraphMaker replaces only the inputs a default names, so the filter is always named: under
  // log_fluorescence the classes that have an SE, under hill_gate_index every class it writes.
  options.push({
    inputName: "filters",
    selectedSource: affinityClass.spec,
    filterType: "equals",
    selectedFilterValues: kdSe ? ["Good", "Partial"] : ["Good", "Partial", "Failed"],
  });
  return options;
});
</script>

<template>
  <TiteseqPage title="Affinity vs Fit Quality" mode="graph">
    <GraphMaker
      v-model="app.model.data.graphStateAffinityVsFit"
      chart-type="scatterplot"
      :data-state-key="app.model.outputs.summaryPfHandle"
      :p-frame="app.model.outputs.summaryPfHandle"
      :default-options="defaultOptions"
      :status-text="{
        noPframe: { title: 'Configure the block in Settings on the Main page and run it.' },
      }"
    />
  </TiteseqPage>
</template>
