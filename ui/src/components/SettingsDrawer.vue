<script setup lang="ts">
import type { DatasetSelection, PlRef } from "@platforma-sdk/model";
import { createDatasetSelection, createPrimaryRef, plRefsEqual } from "@platforma-sdk/model";
import {
  PlAccordionSection,
  PlAlert,
  PlCheckbox,
  PlDatasetSelector,
  PlDropdown,
  PlDropdownRef,
  PlNumberField,
  PlSlideModal,
  PlTooltip,
} from "@platforma-sdk/ui-vue";
import { MIN_CONCENTRATION_POINTS } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed, watch } from "vue";
import { useApp } from "../app";
import { aboveZero, boundedBy, emptyError, wholeNumber } from "../fieldValidation";
import { useFitSettings } from "../composables/useFitSettings";
import { EQUATION_OPTIONS, REPLICATE_ERROR_MODEL_OPTIONS } from "../fitSettings";

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

// The two pairs each field is bounded by the other half of.
const r2GoodAtLeastFailed = boundedBy(
  "≥",
  () => app.model.data.r2ThresholdFailed,
  "R² threshold (Failed)",
);
const r2FailedAtMostGood = boundedBy(
  "≤",
  () => app.model.data.r2ThresholdGood,
  "R² threshold (Good)",
);
const nMinBelowMax = boundedBy("<", () => app.model.data.nMax, "Hill coefficient — max");
const nMaxAboveMin = boundedBy(">", () => app.model.data.nMin, "Hill coefficient — min");

const { setEquation } = useFitSettings();

watch(
  () => app.model.outputs.isRunning,
  (isRunning) => {
    if (isRunning) open.value = false;
  },
);

// Page-level alerts cover spec-based checks only — a sort fraction without a bin column, a
// target antigen without an antigen column. Numeric-field errors show inline on each
// PlNumberField.
const warnings = computed(() => app.model.outputs.validationWarnings ?? []);

const binMode = computed(() => app.model.outputs.binMode === true);

// Every replicate setting rides on the column, so clearing it hides them all. The stored values
// stay: re-picking a column brings the previous combination back.
const hasReplicates = computed(() => app.model.data.replicateColumnRef !== undefined);

// Each equation reads only its own settings. A hidden setting keeps its stored value.
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

// Clear dependent fields in the change handlers, not in a watcher on data. A watcher also fires
// when the webview loads and when another client patches data.

// A target antigen belongs to one antigen column.
function setAntigenColumn(ref: PlRef | undefined): void {
  const prev = app.model.data.antigenColumnRef;
  app.model.data.antigenColumnRef = ref;
  const unchanged = ref !== undefined && prev !== undefined && plRefsEqual(ref, prev);
  if (!unchanged) app.model.data.targetAntigen = undefined;
}

// The model rejects a sort fraction without a bin column.
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
        The Kd,app unit comes from the column label: the text in its trailing parentheses, so
        "Antigen concentration (M)" gives "M". A label without parentheses is used whole.
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
    <PlDropdownRef
      v-model="app.model.data.replicateColumnRef"
      :options="app.model.outputs.replicateOptions"
      label="Sample replicate (optional)"
      clearable
    >
      <template #tooltip>
        Optional. A per-sample column naming the replicate each sample belongs to. Each replicate is
        a titration of its own: it is fitted alone, and the per-clonotype values are then combined.
        The table then shows the combined Kd,app with its standard error and the SD across the
        replicates, and the Replicates page shows each replicate's own Kd,app. Leave empty to fit
        one titration.
      </template>
    </PlDropdownRef>

    <PlAlert v-for="(w, i) in warnings" :key="i" :type="w.severity === 'error' ? 'error' : 'warn'">
      {{ w.message }}
    </PlAlert>

    <PlAccordionSection label="Fit method">
      <PlDropdown
        :model-value="app.model.data.equation"
        :options="EQUATION_OPTIONS"
        label="Equation"
        @update:model-value="(eq) => eq && setEquation(eq)"
      >
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
        :validate="aboveZero"
        clearable
        :error-message="emptyError(app.model.data.plateauFloor)"
      >
        <template #tooltip>
          The plateau that each clonotype must be able to reach at saturation, as a fraction of the
          gate range. The default is 0.667, for a screening library. For a library in which most
          clonotypes bind, use a lower value, for example 0.4.
        </template>
      </PlNumberField>
      <PlCheckbox v-if="isLogFluorescence" v-model="app.model.data.r2GateEnabled">
        Fail in-range fits of poor quality
        <PlTooltip class="info" position="top">
          <template #tooltip>
            When this is on, the block checks each fit that has a Kd,app inside the titrated range.
            A fit gets the rating Failed if its R² is below Minimum R². A fit also gets the rating
            Failed if its standard error is above Maximum Kd,app SE. When this is off, the block
            does not do these checks.
          </template>
        </PlTooltip>
      </PlCheckbox>
      <PlNumberField
        v-if="isLogFluorescence && app.model.data.r2GateEnabled"
        v-model="app.model.data.r2GateMin"
        label="Minimum R²"
        :min-value="0"
        :max-value="1"
        :step="0.05"
        clearable
        :error-message="emptyError(app.model.data.r2GateMin)"
      >
        <template #tooltip>
          The lowest R² that a fit can have. A fit with a lower R², or with no R², gets the rating
          Failed, with the reason low_r2. The block checks only a fit whose Kd,app is inside the
          titrated range. With sample replicates, the block checks the fit of each replicate.
        </template>
      </PlNumberField>
      <PlNumberField
        v-if="isLogFluorescence && app.model.data.r2GateEnabled"
        v-model="app.model.data.kdSeGateMax"
        label="Maximum Kd,app SE (log10, optional)"
        :min-value="0"
        :step="0.1"
        :validate="aboveZero"
        clearable
      >
        <template #tooltip>
          The largest standard error of log10 Kd,app that a fit can have. A fit with a larger error,
          or with no error, gets the rating Failed, with the reason high_kd_se. The block checks
          only a fit whose Kd,app is inside the titrated range. With sample replicates, the block
          checks the fit of each replicate. An error above 0.5 always changes the rating from Good
          to Partial. Leave this field empty to set no Failed limit.
        </template>
      </PlNumberField>
    </PlAccordionSection>

    <!--
      The whole Sample replicates section is hidden while it is being decided what stays. The
      error across replicates is the sample SD over the square root of the replicate count; the
      other model and the two pre-combine filters are not offered. `deriveArgs` pins the error
      model to match, so a block that stored the other one still gets the sample SD. The code
      behind every field here is untouched: uncomment this block and drop that pin to restore it.

    <PlAccordionSection v-if="hasReplicates" label="Sample replicates">
      <PlDropdown
        :model-value="app.model.data.replicateErrorModel"
        :options="REPLICATE_ERROR_MODEL_OPTIONS"
        label="Error across replicates"
        @update:model-value="(m) => m && (app.model.data.replicateErrorModel = m)"
      >
        <template #tooltip>
          How the error of the combined Kd,app is computed. The value itself is the plain mean of
          log10 Kd,app either way, and with one replicate both report that fit's own standard error.
          <br /><br />
          Sample SD / √n: the standard error of the mean, over the replicates that entered. The
          default. <br /><br />
          Population SD / √(n+1): divides by one more than the replicate count and uses the
          population standard deviation, which makes it the smaller of the two. Its own value at one
          replicate is 0, which would rate every single-replicate variant Good, so the fit's own
          error stands in.
        </template>
      </PlDropdown>
      <!--
        The code behind each one is untouched, so uncommenting a field brings it back as it was.

        Each hidden field holds its default, which is what the block now does:
          * Fit error cutoff — unset, so no replicate is dropped for a wide fit error.
          * R² floor — unset, so no replicate is dropped for a poor R². Minimum R², under Fit
            method, already rates such a fit Failed before this floor would see it.
          * Censoring majority — 0.5, so one edge of the titration must simply hold more of the
            replicates than the other for that edge to be reported. The operator asked for this
            to be a fixed "more than half" rather than a number a reader has to choose.

      <PlNumberField
        v-if="isLogFluorescence"
        v-model="app.model.data.replicateFitErrorCutoff"
        label="Fit error cutoff (decades)"
        :min-value="0"
        :step="0.1"
        :validate="aboveZero"
        clearable
      >
        <template #tooltip>
          Optional. Leaves out a replicate whose own fit error is above this, before the mean. Off
          unless you set it. A replicate censored at an edge of the titration is not affected: the
          edge, not the fit, set its value.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.replicateR2Floor"
        label="R² floor"
        :min-value="0"
        :max-value="1"
        :step="0.05"
        clearable
      >
        <template #tooltip>
          Optional. Leaves out a replicate whose fit R² is below this, before the mean, and counts
          it as filtered. Off unless you set it. This differs from Minimum R² under Fit method: that
          gate runs first and rates the replicate's own fit Failed, so a replicate it rejects never
          reaches this floor. Set both to the same value and this one does nothing.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.replicateCensoredMajority"
        label="Censoring majority"
        :min-value="0.5"
        :max-value="1"
        :step="0.05"
        clearable
        :error-message="emptyError(app.model.data.replicateCensoredMajority)"
      >
        <template #tooltip>
          When no replicate measured a Kd,app inside the titrated range, the share of the replicates
          one edge must hold, and exceed the other edge by, for that edge to be reported. Below it
          the clonotype is reported as Failed. This is the block's own rule: no published method
          defines one, and averaging opposite edges would land inside the range and mean nothing.
        </template>
      </PlNumberField>
    </PlAccordionSection>
    -->

    <PlAccordionSection label="Read coverage">
      <PlNumberField
        v-model="app.model.data.minReadsPerConcentration"
        label="Min reads per concentration"
        :min-value="1"
        :step="1"
        :validate="wholeNumber"
        clearable
        :error-message="emptyError(app.model.data.minReadsPerConcentration)"
      >
        <template #tooltip>
          The fewest reads a clonotype needs at one concentration for that point to enter the fit. A
          point below the floor is dropped, and the weighting handles what coverage is left. Both
          equations read it. The default is 3.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.minConcentrationPoints"
        label="Min concentration points"
        :min-value="MIN_CONCENTRATION_POINTS"
        :step="1"
        :validate="wholeNumber"
        clearable
        :error-message="emptyError(app.model.data.minConcentrationPoints)"
      >
        <template #tooltip>
          The minimum number of concentration points that a variant must have for the block to fit
          it. The block counts only the points that remain after the read floor. It does not count
          the 0 M point. A variant with fewer points gets the rating Failed, with the reason
          insufficient_points. The default and the minimum value are both 5. Both equations use this
          setting.
        </template>
      </PlNumberField>
      <!-- Absent reads as off: the field has no migration default. -->
      <PlCheckbox
        :model-value="app.model.data.requireAllConcentrations === true"
        @update:model-value="(on: boolean) => (app.model.data.requireAllConcentrations = on)"
      >
        Require reads at every concentration
        <PlTooltip class="info" position="top">
          <template #tooltip>
            When this is on, a variant with no reads at one or more concentrations gets the rating
            Failed, with the reason incomplete_titration. The 0 M sample counts as a concentration.
            The block does this check before the fit. Both equations use this setting. It is off by
            default. Phillips 2021 fitted only complete titrations.
          </template>
        </PlTooltip>
      </PlCheckbox>
    </PlAccordionSection>

    <PlAccordionSection v-if="!isLogFluorescence" label="Fit quality">
      <PlNumberField
        v-model="app.model.data.r2ThresholdGood"
        label="R² threshold — Good"
        :min-value="0"
        :max-value="1"
        :step="0.05"
        :validate="r2GoodAtLeastFailed"
        clearable
        :error-message="emptyError(app.model.data.r2ThresholdGood)"
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
        :validate="r2FailedAtMostGood"
        clearable
        :error-message="emptyError(app.model.data.r2ThresholdFailed)"
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
        :validate="nMinBelowMax"
        clearable
        :error-message="emptyError(app.model.data.nMin)"
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
        :validate="nMaxAboveMin"
        clearable
        :error-message="emptyError(app.model.data.nMax)"
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
        clearable
        :error-message="emptyError(app.model.data.hookEffectThresholdBin)"
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
        clearable
        :error-message="emptyError(app.model.data.hookEffectThresholdNoBin)"
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
        :validate="wholeNumber"
        clearable
        :error-message="emptyError(app.model.data.hookEffectMinReads)"
      >
        <template #tooltip>
          Skip the hook-effect check when the top two concentration points have fewer reads than
          this. Below the floor, a signal drop is more likely noise than a real hook. Default 20.
        </template>
      </PlNumberField>
    </PlAccordionSection>
  </PlSlideModal>
</template>
