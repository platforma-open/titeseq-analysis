<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { PlAlert, PlTooltip } from "@platforma-sdk/ui-vue";
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

    <div class="conc-unit-hint">
      <PlTooltip position="top">
        <span class="conc-unit-hint__label">
          X-axis is Concentration in attomolar (aM)
          <span class="conc-unit-hint__icon">&#9432;</span>
        </span>
        <template #tooltip>
          Concentrations are stored as integers in attomolar (1 aM = 10<sup>-18</sup> M) so the axis
          is exact across language layers. To read tick values back in molar units:
          <ul>
            <li>10<sup>6</sup> aM = 1 pM (10<sup>-12</sup> M)</li>
            <li>10<sup>9</sup> aM = 1 nM (10<sup>-9</sup> M)</li>
            <li>10<sup>12</sup> aM = 1 µM (10<sup>-6</sup> M)</li>
            <li>10<sup>15</sup> aM = 1 mM (10<sup>-3</sup> M)</li>
          </ul>
          Kd,app values reported on this block are in molar (M).
        </template>
      </PlTooltip>
    </div>

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

<style scoped>
.conc-unit-hint {
  margin: 8px 16px 0;
  font-size: 12px;
  color: var(--txt-03, #666);
}
.conc-unit-hint__label {
  display: inline-flex;
  align-items: center;
  gap: 4px;
  cursor: help;
}
.conc-unit-hint__icon {
  font-size: 14px;
  opacity: 0.7;
}
</style>
