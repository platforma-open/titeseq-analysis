<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { PlAlert } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

const binMode = computed(() => app.model.outputs.binMode === true);
const isEmpty = computed(() => app.model.outputs.isEmpty === true);
const hasResults = computed(() => app.model.outputs.isEmpty === false);

const defaultOptions = computed((): PredefinedGraphOption<"scatterplot">[] | undefined => {
  const pCols = app.model.outputs.titrationCurvesPfCols;
  if (!pCols) return undefined;

  const meanBin = pCols.find((p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/meanBin");
  if (!meanBin) return undefined;

  const concAxis = meanBin.spec.axesSpec.find((a) => a.name === "pl7.app/vdj/concentration");
  const clonotypeAxis = meanBin.spec.axesSpec.find(
    (a) => a.name === "pl7.app/vdj/clonotypeKey" || a.name === "pl7.app/vdj/scClonotypeKey",
  );
  if (!concAxis || !clonotypeAxis) return undefined;

  const affinityClass = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/affinityClass",
  );
  const fittedMeanBin = pCols.find(
    (p: PColumnIdAndSpec) => p.spec.name === "pl7.app/vdj/fittedMeanBin",
  );

  const options: PredefinedGraphOption<"scatterplot">[] = [
    { inputName: "x", selectedSource: concAxis },
    { inputName: "y", selectedSource: meanBin.spec },
    { inputName: "facetBy", selectedSource: clonotypeAxis },
    { inputName: "grouping", selectedSource: clonotypeAxis },
  ];
  if (fittedMeanBin) {
    options.push({ inputName: "additionalCurves", selectedSource: fittedMeanBin.spec });
  }
  if (affinityClass) {
    options.push({
      inputName: "filters",
      selectedSource: affinityClass.spec,
      filterType: "equals",
      selectedFilterValues: ["Good"],
    });
  }
  return options;
});
</script>

<template>
  <TiteseqPage title="Titration Curves">
    <PlAlert v-if="hasResults && !binMode" type="warn">
      No-bin mode: Kd,app values reflect clonotype frequency shifts, not fluorescence. They are not
      comparable to bin-derived results — do not mix in the same Lead Selection ranking.
    </PlAlert>
    <PlAlert v-if="isEmpty" type="warn">
      All clonotypes failed to fit. Loosen thresholds in Inputs or check the Fit log.
    </PlAlert>

    <GraphMaker
      v-model="app.model.data.graphStateTitrationCurves"
      chart-type="scatterplot"
      :data-state-key="app.model.outputs.titrationCurvesPf"
      :p-frame="app.model.outputs.titrationCurvesPf"
      :default-options="defaultOptions"
      :status-text="{
        noPframe: { title: 'Open Inputs (top right) to configure the block and run it.' },
      }"
    />
  </TiteseqPage>
</template>
