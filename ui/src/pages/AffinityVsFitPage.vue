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

  const kdPlot = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/kdPlotPosition");
  const hillPlot = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/hillPlotPosition",
  );
  const hillCoef = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/hillCoefficient",
  );
  const affinityClass = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/affinityClass",
  );
  const failureReason = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/fitFailureReason",
  );
  const kdSe = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/kdLog10Se");
  // y is the measure each equation rates fit quality on: SE of log10 Kd,app, or Hill coefficient.
  const ySource = app.model.data.equation === "log_fluorescence" && kdSe ? kdSe : hillPlot;
  if (!kdPlot || !ySource || !failureReason) return undefined;

  // A failed or pinned log-fluorescence fit has no SE. GraphMaker replaces only the inputs a
  // default names, so the filter is named here to replace the saved Failed-only filter.
  if (ySource === kdSe) {
    const options: PredefinedGraphOption<"scatterplot">[] = [
      { inputName: "x", selectedSource: kdPlot.spec },
      { inputName: "y", selectedSource: kdSe.spec },
    ];
    if (affinityClass) {
      options.push(
        { inputName: "grouping", selectedSource: affinityClass.spec },
        {
          inputName: "filters",
          selectedSource: affinityClass.spec,
          filterType: "equals",
          selectedFilterValues: ["Good", "Partial"],
        },
      );
    }
    return options;
  }

  const options: PredefinedGraphOption<"scatterplot">[] = [
    { inputName: "x", selectedSource: kdPlot.spec },
    { inputName: "y", selectedSource: ySource.spec },
    { inputName: "grouping", selectedSource: failureReason.spec },
  ];
  if (ySource === hillPlot && hillCoef) {
    options.push({
      inputName: "filters",
      selectedSource: hillCoef.spec,
      filterType: "range",
      selectedFilterRange: { min: 0 },
    });
  }
  if (affinityClass) {
    options.push({
      inputName: "filters",
      selectedSource: affinityClass.spec,
      filterType: "equals",
      selectedFilterValues: ["Failed"],
    });
  }
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
        noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
      }"
    />
  </TiteseqPage>
</template>
