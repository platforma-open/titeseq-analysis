<script setup lang="ts">
import type { DatasetSelection, PlRef } from "@platforma-sdk/model";
import { createDatasetSelection, createPrimaryRef, plRefsEqual } from "@platforma-sdk/model";
import {
  PlAccordionSection,
  PlAlert,
  PlBtnGhost,
  PlElementList,
  PlMaskIcon24,
  PlDatasetSelector,
  PlDropdown,
  PlDropdownRef,
  PlNumberField,
  PlSlideModal,
  PlTooltip,
} from "@platforma-sdk/ui-vue";
import {
  MIN_CONCENTRATION_POINTS,
  refKey,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed, ref, watch } from "vue";
import { useApp } from "../app";
import { closedSettingsWarnings } from "../closedSettingsWarnings";
import { aboveZero, boundedBy, emptyError, wholeNumber } from "../fieldValidation";
import { useFitSettings } from "../composables/useFitSettings";
import { useNoun } from "../composables/useNoun";
import { EQUATION_OPTIONS } from "../fitSettings";

const open = defineModel<boolean>({ required: true });

const app = useApp();
const { noun, Noun } = useNoun();

// The selector picks an abundance, or an abundance narrowed by one of its subset columns.
// Stored as two refs so blocks saved before filters existed keep their `abundanceRef`.
const abundanceSelection = computed<DatasetSelection | undefined>({
  get: () => {
    const { abundanceRef, filterRef } = app.model.data;
    if (abundanceRef === undefined) return undefined;
    return createDatasetSelection(createPrimaryRef(abundanceRef, filterRef));
  },
  // The variant axis's alphabet is snapshotted with the pick: the args cannot read the result
  // pool. Unknown while the output is not loaded.
  set: (selection) => {
    const column = selection?.primary.column;
    const nucleotide = app.model.outputs.nucleotideAbundances;
    app.model.data.abundanceRef = column;
    app.model.data.filterRef = selection?.primary.filter;
    app.model.data.variantAxisIsNucleotide =
      column === undefined || nucleotide === undefined
        ? undefined
        : nucleotide.includes(refKey(column));
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

// Advice about the settings. The user can close each one.
const settingsWarnings = computed(() =>
  (app.model.outputs.settingsWarnings ?? []).filter((m) => !closedSettingsWarnings.value.has(m)),
);
const closeSettingsWarning = (message: string) => {
  closedSettingsWarnings.value = new Set([...closedSettingsWarnings.value, message]);
};

const binMode = computed(() => app.model.outputs.binMode === true);

// Each equation reads only its own settings. A hidden setting keeps its stored value.
const isLogFluorescence = computed(() => app.model.data.equation === "log_fluorescence");

// The log fluorescence fit cannot run without bins. Run stays disabled until one is picked, so the
// field says why.
const binMissing = computed(
  () => isLogFluorescence.value && app.model.data.binColumnRef === undefined,
);

// The bin-dependent fields: shown when a bin is picked, or when the fit needs one.
const readsBins = computed(
  () => isLogFluorescence.value || app.model.data.binColumnRef !== undefined,
);

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

const binFacts = (ref: PlRef | undefined) =>
  ref === undefined ? undefined : app.model.outputs.binFactsByRef?.[refKey(ref)];

/** The picked bin column's values in their default order: numeric, else alphabetical. */
const binValues = (ref: PlRef | undefined): string[] | undefined => binFacts(ref)?.values;

// The model rejects a sort fraction without a bin column. The order and the column's type are
// snapshotted in the same gesture as the pick: an order from another column means nothing for
// this one, and the Run gate reads the type from data.
function setBinColumn(ref: PlRef | undefined): void {
  const facts = binFacts(ref);
  app.model.data.binColumnRef = ref;
  app.model.data.binOrder = facts?.values;
  app.model.data.binColumnIsText = facts?.isText;
  if (ref === undefined) app.model.data.sortFractionColumnRef = undefined;
}

// Bin values the column holds and the stored order does not rank. The run refuses them.
const unrankedBins = computed(() => app.model.outputs.unrankedBins ?? []);

// A label column's default order is alphabetical, which is rarely the gate order, so its section
// starts open. A number column's default is its own numbering, so its section starts closed.
// View state: it lives in the component and resets with it.
const binOrderOpenChoice = ref<boolean>();
const binOrderOpen = computed({
  get: () =>
    binOrderOpenChoice.value ??
    (app.model.data.binColumnIsText === true || unrankedBins.value.length > 0),
  set: (value: boolean) => {
    binOrderOpenChoice.value = value;
  },
});

// A user gesture: the new values go last, as the highest gates, and the user drags them into place.
function addUnrankedBins(): void {
  const order = app.model.data.binOrder;
  if (order === undefined) return;
  app.model.data.binOrder = [...order, ...unrankedBins.value];
}

const binOrderModified = computed(() => {
  const current = app.model.data.binOrder;
  const base = binValues(app.model.data.binColumnRef);
  if (current === undefined || base === undefined) return false;
  return current.length !== base.length || current.some((v, i) => v !== base[i]);
});

// A user gesture, so writing data from an output here is not the hairpin.
function resetBinOrder(): void {
  app.model.data.binOrder = binValues(app.model.data.binColumnRef);
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
        Per-sample, per-{{ noun }} integer read counts. Selecting this column also anchors the block
        to that upstream dataset. Normalized per-sample inside the block. On a nucleotide variant
        dataset whose amino-acid linker is present, each nucleotide variant is fitted on its own and
        the fits are combined per amino-acid variant.
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
    <!-- Above FACS bin: the equation decides whether a bin is required and which bin fields show. -->
    <PlDropdown
      :model-value="app.model.data.equation"
      :options="EQUATION_OPTIONS"
      label="Fit method equation"
      @update:model-value="(eq) => eq && setEquation(eq)"
    >
      <template #tooltip>
        Log fluorescence fit: fits a Hill curve to the mean log10 fluorescence of each {{ noun }}.
        The Hill exponent is fixed at 1. New blocks use this fit. <br /><br />
        Mean bin fit: fits a Hill curve to the read-weighted mean bin of each {{ noun }}. The Hill
        exponent is free. Blocks created before the equation choice existed keep this fit.
        <br /><br />
        Each equation rates fit quality in its own way. "Good" means that the equation found a
        Kd,app inside the titrated range. It does not give an error bound.
      </template>
    </PlDropdown>
    <PlDropdownRef
      :model-value="app.model.data.binColumnRef"
      :options="app.model.outputs.binOptions"
      @update:model-value="setBinColumn"
      :label="isLogFluorescence ? 'FACS bin' : 'FACS bin (optional)'"
      :required="isLogFluorescence"
      :error="binMissing ? 'The Log fluorescence fit needs a FACS bin column.' : undefined"
      clearable
    >
      <template #tooltip>
        A per-sample value that names the FACS bin: a number or a label. Set the bins' order below.
        The Log fluorescence fit needs it. The Mean bin fit can run without it, in no-bin mode.
        Kd,app from no-bin mode is not comparable to Kd,app from bins.
      </template>
    </PlDropdownRef>
    <!-- The gates from low to high fluorescence. Position is the rank, so the values can be any
         labels. No value can be removed: a gate left out would bias every mean bin. A block with
         no stored order reads the column's integers until the user sets one. -->
    <template v-if="app.model.data.binColumnRef !== undefined">
      <PlAccordionSection
        v-if="app.model.data.binOrder !== undefined"
        v-model="binOrderOpen"
        label="FACS bin order"
      >
        <div class="titeseq-bin-order-label">
          Lowest fluorescence at the top
          <PlTooltip class="info">
            <template #tooltip>
              Drag the bins into order: the lowest fluorescence at the top, the highest at the
              bottom. The position sets each bin's rank, so the bin values do not need to be
              numbers. Every gate takes part: leaving one out would drop its cells and bias every
              mean bin. An unsorted sample has an empty bin value.
            </template>
          </PlTooltip>
        </div>
        <PlElementList v-model:items="app.model.data.binOrder" disable-removing>
          <template #item-title="{ item }">{{ item }}</template>
        </PlElementList>
        <PlAlert v-if="unrankedBins.length > 0" type="error">
          The bin column has values that the order does not rank: {{ unrankedBins.join(", ") }}. The
          run fails until every bin is ranked.
          <PlBtnGhost @click="addUnrankedBins">Add to the order</PlBtnGhost>
        </PlAlert>
        <PlBtnGhost v-if="binOrderModified" @click="resetBinOrder">
          Reset to default
          <template #append>
            <PlMaskIcon24 name="reverse" />
          </template>
        </PlBtnGhost>
      </PlAccordionSection>
      <template v-else>
        <PlAlert v-if="app.model.data.binColumnIsText === true" type="error">
          The bin column holds labels, not numbers, so the block cannot run until the bin order is
          set.
          <template v-if="binValues(app.model.data.binColumnRef) === undefined">
            The column's values are not loaded yet.
          </template>
        </PlAlert>
        <PlBtnGhost
          v-if="binValues(app.model.data.binColumnRef) !== undefined"
          @click="resetBinOrder"
        >
          Set the bin order
        </PlBtnGhost>
      </template>
    </template>
    <!-- Shown whenever the fit reads bins, so the fields a required bin brings are visible before
         the bin is picked. -->
    <PlDropdownRef
      v-if="readsBins"
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
      v-if="isLogFluorescence"
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
        a titration of its own: it is fitted alone, and the per-{{ noun }} values are then combined.
        The table then shows the combined Kd,app with its standard error and the SD across the
        replicates, and the Replicates page shows each replicate's own Kd,app. Leave empty to fit
        one titration.
      </template>
    </PlDropdownRef>

    <PlAlert v-for="(w, i) in warnings" :key="i" :type="w.severity === 'error' ? 'error' : 'warn'">
      {{ w.message }}
    </PlAlert>
    <PlAlert
      v-for="m in settingsWarnings"
      :key="m"
      type="warn"
      closeable
      @update:model-value="(shown: boolean) => shown || closeSettingsWarning(m)"
    >
      {{ m }}
    </PlAlert>

    <!-- Every field here belongs to the log fluorescence fit; the mean bin fit has its own section. -->
    <PlAccordionSection v-if="isLogFluorescence" label="Fit method">
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
          The plateau that each {{ noun }} must be able to reach at saturation, as a fraction of the
          gate range. The default is 0.667, for a screening library. For a library in which most
          {{ noun }}s bind, use a lower value, for example 0.4.
        </template>
      </PlNumberField>
      <PlNumberField
        v-if="isLogFluorescence"
        v-model="app.model.data.r2GateMin"
        label="Minimum R² (optional)"
        :min-value="0"
        :max-value="1"
        :step="0.05"
        clearable
      >
        <template #tooltip>
          The lowest R² that a fit can have. A fit with a lower R², or with no R², gets the rating
          Failed, with the reason low_r2. The block checks only a fit whose Kd,app is inside the
          titrated range. With sample replicates, the block checks the fit of each replicate. The
          default is 0.8. Leave this field empty to set no R² limit.
        </template>
      </PlNumberField>
      <PlNumberField
        v-if="isLogFluorescence"
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
          to Partial. The default is 1. Leave this field empty to set no Failed limit.
        </template>
      </PlNumberField>
    </PlAccordionSection>

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
          The fewest reads a {{ noun }} needs at one concentration for that point to enter the fit.
          A point below the floor is dropped, and the weighting handles what coverage is left. Both
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
      <!-- "Require reads at every concentration" is withheld: the run always fits a titration
        with a gap. -->
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
          {{ Noun }}s with weighted R² at or above this threshold and Hill coefficient in range
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
          {{ Noun }}s with weighted R² below this threshold classify as Failed regardless of n.
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
          the top concentration's {{ noun }} frequency drops by more than this.
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

<style scoped>
.titeseq-bin-order-label {
  display: flex;
  margin-bottom: -15px;
}
</style>
