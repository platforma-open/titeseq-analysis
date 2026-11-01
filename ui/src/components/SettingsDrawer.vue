<script setup lang="ts">
import type { DatasetSelection, PlRef } from "@platforma-sdk/model";
import { createDatasetSelection, createPrimaryRef, plRefsEqual } from "@platforma-sdk/model";
import {
  PlAccordionSection,
  PlAlert,
  PlDatasetSelector,
  PlDropdown,
  PlDropdownRef,
  PlNumberField,
  PlSlideModal,
} from "@platforma-sdk/ui-vue";
import { computed, watch } from "vue";
import { useApp } from "../app";
import { useFieldValidation } from "../composables/useFieldValidation";
import { useFitSettings } from "../composables/useFitSettings";
import { EQUATION_OPTIONS } from "../fitSettings";

const open = defineModel<boolean>({ required: true });

const app = useApp();

// The selector picks an abundance, or an abundance narrowed by one of its subset columns.
// Stored as two refs so blocks saved before filters existed keep their `abundanceRef`.
const abundanceSelection = computed<DatasetSelection | undefined>({
  get: () => {
    const { abundanceRef, filterRef } = app.model.data;
    if (abundanceRef === undefined) return undefined;
    return createDatasetSelection(createPrimaryRef(abundanceRef, filterRef));
  },
  set: (selection) => {
    app.model.data.abundanceRef = selection?.primary.column;
    app.model.data.filterRef = selection?.primary.filter;
  },
});

const {
  plateauFloorError,
  minReadsError,
  minConcPointsError,
  r2GoodError,
  r2FailedError,
  nMinError,
  nMaxError,
  hookThresholdBinError,
  hookThresholdNoBinError,
  hookMinReadsError,
} = useFieldValidation();

const { presetOptions, selectedPreset } = useFitSettings();

watch(
  () => app.model.outputs.isRunning,
  (isRunning) => {
    if (isRunning) open.value = false;
  },
);

// Page-level alerts cover spec-based checks only — concentration column
// label format, sort fraction without bin column. Numeric-field bounds show
// inline on each PlNumberField.
const warnings = computed(() => app.model.outputs.validationWarnings ?? []);

const binMode = computed(() => app.model.outputs.binMode === true);

// The log-fluorescence fit never reads the v2.4.0 settings, and the mean bin fit never reads the
// plateau floor. Hiding is display-only: the stored values stay as they are.
const isLogFluorescence = computed(() => app.model.data.equation === "log_fluorescence");

const targetAntigenOptions = computed(() =>
  (app.model.outputs.targetAntigenValues ?? []).map((v) => ({ value: v, label: v })),
);

// Track the selected antigen column's human label so the value-picker label
// reads as "Target <column>" (e.g. "Target Sample" when the user chose a
// metadata column labelled "Sample") rather than the fixed "Target antigen".
const selectedAntigenColumnLabel = computed(() => {
  const ref = app.model.data.antigenColumnRef;
  if (!ref) return undefined;
  const options = app.model.outputs.antigenOptions ?? [];
  const match = options.find((o) => o.ref.blockId === ref.blockId && o.ref.name === ref.name);
  return match?.label;
});

const targetAntigenLabel = computed(() =>
  selectedAntigenColumnLabel.value
    ? `Target ${selectedAntigenColumnLabel.value}`
    : "Target antigen",
);

// Dependent fields are cleared in the pickers' own change handlers, never by a watcher on data:
// a watcher also fires when the webview loads or another client patches data, and wiped
// targetAntigen that way.

// A target antigen belongs to one antigen column.
function setAntigenColumn(ref: PlRef | undefined): void {
  const prev = app.model.data.antigenColumnRef;
  app.model.data.antigenColumnRef = ref;
  const unchanged = ref !== undefined && prev !== undefined && plRefsEqual(ref, prev);
  if (!unchanged) app.model.data.targetAntigen = undefined;
}

// Sort-fraction correction applies in bin mode only. The model rejects a sort fraction without a
// bin column, so clearing the bin column clears the sort fraction too.
function setBinColumn(ref: PlRef | undefined): void {
  app.model.data.binColumnRef = ref;
  if (ref === undefined) app.model.data.sortFractionColumnRef = undefined;
}
</script>

<template>
  <PlSlideModal v-model="open" close-on-outside-click shadow>
    <template #title>Settings</template>

    <PlDatasetSelector
      v-model="abundanceSelection"
      :options="app.model.outputs.abundanceOptions"
      label="Read count column"
      required
    >
      <template #tooltip>
        Per-sample, per-clonotype integer read counts from MiXCR. Selecting this column also anchors
        the block to that upstream dataset. Normalized per-sample inside the block.
      </template>
    </PlDatasetSelector>
    <PlDropdownRef
      v-model="app.model.data.concentrationColumnRef"
      :options="app.model.outputs.concentrationOptions"
      label="Antigen concentration"
      required
    >
      <template #tooltip>
        Per-sample numeric column giving the antigen concentration at which each sample was stained.
        The column label becomes the Kd,app unit — prefer a bare unit (e.g. "nM", "µM") over a
        phrase.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      :model-value="app.model.data.binColumnRef"
      :options="app.model.outputs.binOptions"
      @update:model-value="setBinColumn"
      label="FACS bin"
      clearable
    >
      <template #tooltip>
        Per-sample positive integer identifying the FACS bin. Leave empty to run in no-bin mode;
        resulting Kd,app values are not comparable to bin-derived values.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      v-if="app.model.data.binColumnRef !== undefined"
      v-model="app.model.data.sortFractionColumnRef"
      :options="app.model.outputs.sortFractionOptions"
      label="FACS sort fraction (optional)"
      clearable
    >
      <template #tooltip>
        Optional. Per-sample numerical metadata column giving the fraction of cells sorted into that
        sample's (concentration, bin) — C_bc / C_c in Adams, Mora, Walczak, Kinney 2016. Values must
        sum to 1 per concentration. When supplied, Mean Bin is corrected for FACS sort yield.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      v-if="isLogFluorescence && app.model.data.binColumnRef !== undefined"
      v-model="app.model.data.binFluorescenceColumnRef"
      :options="app.model.outputs.binFluorescenceOptions"
      label="Gate fluorescence (optional)"
      clearable
    >
      <template #tooltip>
        Optional. A per-sample number column with the measured mean log10 fluorescence of each gate.
        If you do not select a column, the fit assumes equal distances between the gates. This
        assumption is less accurate for a real sorter.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      :model-value="app.model.data.antigenColumnRef"
      :options="app.model.outputs.antigenOptions"
      @update:model-value="setAntigenColumn"
      label="Antigen label"
      clearable
    >
      <template #tooltip>
        Optional per-sample antigen name. Select when a dataset contains samples stained against
        multiple antigens and you want to analyse one of them.
      </template>
    </PlDropdownRef>
    <PlDropdown
      v-if="app.model.data.antigenColumnRef"
      v-model="app.model.data.targetAntigen"
      :options="targetAntigenOptions"
      :label="targetAntigenLabel"
      :error-status="!app.model.data.targetAntigen"
      required
    >
      <template #tooltip>
        Which antigen to analyse. Required when an antigen column is set.
      </template>
    </PlDropdown>

    <PlAlert v-for="(w, i) in warnings" :key="i" :type="w.severity === 'error' ? 'error' : 'warn'">
      {{ w.message }}
    </PlAlert>

    <PlAccordionSection label="Fit method">
      <PlDropdown
        v-model="selectedPreset"
        :options="presetOptions"
        label="Published study preset"
        placeholder="None"
        clearable
      >
        <template #tooltip>
          Sets the fit settings to the values that a published study used. You can change each
          setting after you select a preset. The selector then shows "(modified)". To use the
          published values again, clear the preset and select it again.
        </template>
      </PlDropdown>
      <PlDropdown v-model="app.model.data.equation" :options="EQUATION_OPTIONS" label="Equation">
        <template #tooltip>
          Log fluorescence fit: fits a Hill curve to the mean log10 fluorescence of each clonotype.
          The Hill exponent is fixed at 1. New blocks use this fit.
          <br /><br />
          Mean bin fit: fits a Hill curve to the read-weighted mean bin of each clonotype. The Hill
          exponent is free. Blocks created before the equation choice existed keep this fit.
          <br /><br />
          Each equation rates fit quality in its own way. "Good" means that the equation found a
          Kd,app inside the titrated range. It does not give an error bound.
        </template>
      </PlDropdown>
      <PlNumberField
        v-if="isLogFluorescence"
        v-model="app.model.data.plateauFloor"
        label="Plateau floor"
        :min-value="0"
        :max-value="1"
        :step="0.05"
        :error-message="plateauFloorError"
        required
      >
        <template #tooltip>
          The plateau that each clonotype must be able to reach at saturation, as a fraction of the
          gate range. The default is 2/3, for a screening library. For a library in which most
          clonotypes bind, use a lower value, for example 0.4.
        </template>
      </PlNumberField>
    </PlAccordionSection>

    <PlAccordionSection v-if="!isLogFluorescence" label="Read coverage">
      <PlNumberField
        v-model="app.model.data.minReadsPerConcentration"
        label="Min reads per concentration"
        :min-value="1"
        :step="1"
        :error-message="minReadsError"
        required
      >
        <template #tooltip>
          Floor applied per clonotype per concentration. Points below the floor are excluded;
          weighting handles remaining coverage variation. Default 3.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.minConcentrationPoints"
        label="Min concentration points"
        :min-value="3"
        :step="1"
        :error-message="minConcPointsError"
        required
      >
        <template #tooltip>
          Minimum number of concentration points (after floor filtering) to attempt a fit. Default 5
          — the practical minimum for a 4-parameter Hill fit.
        </template>
      </PlNumberField>
    </PlAccordionSection>

    <PlAccordionSection v-if="!isLogFluorescence" label="Fit quality">
      <PlNumberField
        v-model="app.model.data.r2ThresholdGood"
        label="R² threshold — Good"
        :min-value="0"
        :max-value="1"
        :step="0.05"
        :error-message="r2GoodError"
        required
      >
        <template #tooltip>
          Clonotypes with weighted R² at or above this threshold and Hill coefficient in range
          classify as Good. Default 0.8.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.r2ThresholdFailed"
        label="R² threshold — Failed"
        :min-value="0"
        :max-value="1"
        :step="0.05"
        :error-message="r2FailedError"
        required
      >
        <template #tooltip>
          Clonotypes with weighted R² below this threshold classify as Failed regardless of n.
          Default 0.5.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.nMin"
        label="Hill coefficient — min"
        :min-value="0"
        :step="0.1"
        :error-message="nMinError"
        required
      >
        <template #tooltip>
          Lower bound on n for Good classification. Values below this downgrade the class by one.
          Default 0.5.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.nMax"
        label="Hill coefficient — max"
        :min-value="0"
        :step="0.1"
        :error-message="nMaxError"
        required
      >
        <template #tooltip>
          For multimeric antigens (trimeric Spike, dimeric receptors), genuine cooperativity can
          push n into 2–4; raise the max accordingly.
        </template>
      </PlNumberField>
    </PlAccordionSection>

    <PlAccordionSection v-if="!isLogFluorescence" label="Hook effect">
      <PlNumberField
        v-if="binMode"
        v-model="app.model.data.hookEffectThresholdBin"
        label="Signal drop threshold (bin mode)"
        :min-value="0"
        :step="0.05"
        :error-message="hookThresholdBinError"
        required
      >
        <template #tooltip>
          Flags non-monotonic signals (potentially genuine tight binders showing a hook effect) when
          the top concentration drops by more than this.
        </template>
      </PlNumberField>
      <PlNumberField
        v-else
        v-model="app.model.data.hookEffectThresholdNoBin"
        label="Signal drop threshold (frequency mode)"
        :min-value="0"
        :step="0.005"
        :error-message="hookThresholdNoBinError"
        required
      >
        <template #tooltip>
          Flags non-monotonic signals (potentially genuine tight binders showing a hook effect) when
          the top concentration's clonotype frequency drops by more than this.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.hookEffectMinReads"
        label="Min reads for hook check"
        :min-value="0"
        :step="1"
        :error-message="hookMinReadsError"
        required
      >
        <template #tooltip>
          Skip the hook-effect check when the top two concentration points have fewer reads than
          this. Below the floor, a signal drop is more likely noise than a real hook. Default 20.
        </template>
      </PlNumberField>
    </PlAccordionSection>
  </PlSlideModal>
</template>
