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
  // Each equation's fit quality on y: the SE of log10 Kd,app for log fluorescence, the Hill
  // coefficient for the mean bin fit. Rows with no SE (failed or pinned fits) do not plot.
  // GraphMaker re-applies these defaults when they change, so switching equation moves y.
  const ySource = app.model.data.equation === "log_fluorescence" && kdSe ? kdSe : hillPlot;
  if (!kdPlot || !ySource || !failureReason) return undefined;

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
