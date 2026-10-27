<script setup lang="ts">
import {
  PlAccordionSection,
  PlAlert,
  PlCheckbox,
  PlDropdown,
  PlDropdownRef,
  PlNumberField,
  PlSectionSeparator,
  PlSlideModal,
} from "@platforma-sdk/ui-vue";
import { plRefsEqual } from "@platforma-sdk/model";
import { computed, watch } from "vue";
import { useApp } from "../app";
import { useFieldValidation } from "../composables/useFieldValidation";

const open = defineModel<boolean>({ required: true });

const app = useApp();

const {
  amplitudePriorError,
  minReadsError,
  minConcPointsError,
  r2GoodError,
  r2FailedError,
  nMinError,
  nMaxError,
  hookThresholdBinError,
  hookThresholdNoBinError,
  hookMinReadsError,
  binWidthError,
  plateauError,
  kdSdGoodError,
  baselineError,
} = useFieldValidation();

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

// Settings follow the unsaved bin choice; `outputs.binMode` is the mode of the results on screen.
const binMode = computed(() => app.model.data.binColumnRef !== undefined);

// The model the block routed to, from the model layer, so the drawer and the fit cannot disagree.
// Undefined only before the first output lands; the ranking fit is the safe thing to show then.
const resolvedModel = computed(() => app.model.outputs.resolvedFitModel ?? "hill_free_n");
const adamsModel = computed(() => resolvedModel.value === "poisson_log");

const modelDescription = computed(() =>
  adamsModel.value
    ? "Adams 2016 measurement model: Poisson likelihood on reads, log fluorescence, Hill n = 1."
    : "Ranking fit: linear Hill on mean bin. Orders candidates; does not measure them.",
);

const fitModelOverrideOptions = [
  { value: undefined, label: "Automatic (from the columns bound)" },
  { value: "poisson_log" as const, label: "Adams 2016 (log fluorescence, n = 1)" },
  { value: "hill_free_n" as const, label: "Linear Hill" },
];

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

// An external data patch replaces the ref with an equal copy; only a different column clears the target.
watch(
  () => app.model.data.antigenColumnRef,
  (next, prev) => {
    if (next !== undefined && prev !== undefined && plRefsEqual(next, prev)) return;
    app.model.data.targetAntigen = undefined;
  },
);

// Sort fractions and cell counts only apply in bin mode. Clearing the bin column
// would leave orphaned refs that the model validator would then reject, so drop
// them here in lock-step with the bin column.
watch(
  () => app.model.data.binColumnRef,
  (binRef) => {
    if (binRef === undefined) {
      app.model.data.sortFractionColumnRef = undefined;
      app.model.data.cellCountColumnRef = undefined;
      app.model.data.binFluorescenceColumnRef = undefined;
      app.model.data.gateWidthColumnRef = undefined;
    }
  },
);
</script>

<template>
  <PlSlideModal v-model="open" close-on-outside-click shadow>
    <template #title>Settings</template>

    <PlDropdownRef
      v-model="app.model.data.abundanceRef"
      :options="app.model.outputs.abundanceOptions"
      label="Read count column"
      required
    >
      <template #tooltip>
        Per-sample, per-clonotype integer read counts from MiXCR. Selecting this column also anchors
        the block to that upstream dataset. Normalized per-sample inside the block.
      </template>
    </PlDropdownRef>
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
      v-model="app.model.data.binColumnRef"
      :options="app.model.outputs.binOptions"
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
      v-model="app.model.data.cellCountColumnRef"
      :options="app.model.outputs.cellCountOptions"
      label="Sorted cell counts"
      clearable
    >
      <template #tooltip>
        Per-sample number of cells FACS sorted into that sample's (concentration, bin) — C_bc in
        Adams, Mora, Walczak, Kinney 2016. Weights each bin by the cells it holds and caps a bin's
        reads at its cell count. Optional: supplying it switches the block from ranking candidates
        to measuring them.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      v-if="
        app.model.data.binColumnRef !== undefined && app.model.data.cellCountColumnRef === undefined
      "
      v-model="app.model.data.sortFractionColumnRef"
      :options="app.model.outputs.sortFractionOptions"
      label="FACS sort fraction (optional)"
      clearable
    >
      <template #tooltip>
        Use when only fractions are known. Per-sample fraction of cells sorted into that sample's
        (concentration, bin) — C_bc / C_c in Adams, Mora, Walczak, Kinney 2016. Values must sum to 1
        per concentration. Sorted cell counts are preferred: fractions cannot cap a bin's reads.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      v-model="app.model.data.antigenColumnRef"
      :options="app.model.outputs.antigenOptions"
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

    <PlAccordionSection v-if="binMode" label="FACS refinements (optional)">
      <PlDropdownRef
        v-model="app.model.data.binFluorescenceColumnRef"
        :options="app.model.outputs.binFluorescenceOptions"
        label="Measured bin fluorescence"
        clearable
      >
        <template #tooltip>
          Per-sample mean log10 fluorescence of the cells sorted into that gate, from the sort
          report. Bound, each bin's fluorescence is read per concentration instead of from one fixed
          gate spacing. A gate's real mean moves with concentration, because at low concentration
          the cells in the high gates sit at the bottom of the gate, and a fixed value makes every
          midpoint read tighter than it is.
        </template>
      </PlDropdownRef>
      <PlDropdownRef
        v-model="app.model.data.gateWidthColumnRef"
        :options="app.model.outputs.gateWidthOptions"
        label="FACS gate width (decades)"
        clearable
      >
        <template #tooltip>
          Per-sample width of one FACS gate in log10 fluorescence. Bound, it supplies the bin width
          below. Every sort report carries it.
        </template>
      </PlDropdownRef>
    </PlAccordionSection>

    <PlAccordionSection v-if="binMode" label="Curve model">
      <PlAlert type="info">{{ modelDescription }}</PlAlert>
      <PlDropdown
        v-model="app.model.data.fitModelOverride"
        :options="fitModelOverrideOptions"
        label="Model"
      >
        <template #tooltip>
          Automatic picks the measurement model when the cells sorted into each sample are bound and
          the ranking fit otherwise. Feeding uncorrected bins to the measurement model is measurably
          worse than ranking them, which is why more data has to switch the model rather than only
          refine it. Override to reproduce an older project's numbers by hand.
        </template>
      </PlDropdown>
      <template v-if="adamsModel">
        <PlNumberField
          v-if="app.model.data.gateWidthColumnRef === undefined"
          v-model="app.model.data.binWidthDecades"
          label="Bin width (decades of fluorescence)"
          :min-value="0"
          :step="0.05"
          :error-message="binWidthError"
          required
        >
          <template #tooltip>
            Width of one FACS gate in log10 fluorescence — the spacing between bin gates. Adams 2016
            used 0.61–0.71. Default 0.65.
          </template>
        </PlNumberField>
        <PlNumberField
          v-model="app.model.data.minPlateauFraction"
          label="Minimum plateau (fraction of bin range)"
          :min-value="0"
          :max-value="1"
          :step="0.05"
          :error-message="plateauError"
          required
        >
          <template #tooltip>
            Every variant is assumed able to reach at least this fraction of the bin range at
            saturation, so a flat curve reads as Kd above the titration rather than a tight binder.
            A tight binder that displays poorly is then reported weaker than it is. Adams 2016 used
            0.63–0.70. Default 0.65.
          </template>
        </PlNumberField>
        <PlNumberField
          v-model="app.model.data.kdSdGood"
          label="Good class: max Kd,app uncertainty (decades)"
          :min-value="0"
          :step="0.05"
          :error-message="kdSdGoodError"
          required
        >
          <template #tooltip>
            An in-range Kd,app is Good when its posterior standard deviation is at or below this,
            Partial above it. Kd,app above the titration range is Failed; below it is Good. Default
            0.5.
          </template>
        </PlNumberField>
        <PlNumberField
          v-model="app.model.data.baselineFluorescence"
          label="0 M baseline fluorescence (optional)"
          :min-value="0"
          :step="0.1"
          :error-message="baselineError"
          clearable
        >
          <template #tooltip>
            Mean fluorescence of cells at 0 M antigen, in units of the lowest bin's fluorescence.
            Leave empty to compute it from the 0 M samples, which are then required.
          </template>
        </PlNumberField>
      </template>
    </PlAccordionSection>

    <PlAccordionSection label="Advanced fit">
      <PlCheckbox v-model="app.model.data.pinHillExponent">
        Fix the Hill exponent at 1
        <template #tooltip>
          The Adams 2016 model fixes it at 1. Fitting it instead costs rank agreement — 0.909 down
          to 0.816 on Adams 2016 replicate 1 — and is what makes the hook-effect check necessary.
          Clear it only to inspect the exponent itself.
        </template>
      </PlCheckbox>
      <PlCheckbox v-model="app.model.data.hookDetectorEnabled">
        Discard non-monotonic (hook effect) curves
        <template #tooltip>
          Off by default. With the exponent fixed at 1 the curves it catches fail to converge on
          their own, and leaving it on discarded 44, 99 and 73 of the 250 tightest binders on the
          three Adams 2016 replicates while buying no rank agreement.
        </template>
      </PlCheckbox>
      <PlNumberField
        v-model="app.model.data.amplitudePriorQ"
        label="Amplitude floor (fraction of bin range)"
        :min-value="0"
        :max-value="1"
        :step="0.05"
        :error-message="amplitudePriorError"
        clearable
      >
        <template #tooltip>
          Empty by default. Raising it makes every variant assumed able to reach this fraction of
          the bin range at saturation. It deletes rows rather than reordering them: a genuinely
          tight binder that displays poorly has a real plateau that is simply low, and raising the
          floor cut 21, 13 and 3 of the 250 tightest binders on the three Adams 2016 replicates.
        </template>
      </PlNumberField>
    </PlAccordionSection>

    <PlAccordionSection label="Read coverage">
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

    <PlAccordionSection v-if="!adamsModel" label="Fit quality">
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

    <PlAccordionSection
      v-if="!adamsModel && app.model.data.hookDetectorEnabled"
      label="Hook effect"
    >
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
