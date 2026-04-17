<script setup lang="ts">
import type { PlRef } from "@platforma-sdk/model";
import {
  PlAccordionSection,
  PlAlert,
  PlDropdownRef,
  PlNumberField,
  PlSectionSeparator,
  PlSlideModal,
  PlTextField,
} from "@platforma-sdk/ui-vue";
import { computed, watch } from "vue";
import { useApp } from "../app";

const open = defineModel<boolean>({ required: true });

const app = useApp();

watch(
  () => app.model.outputs.isRunning,
  (isRunning) => {
    if (isRunning) open.value = false;
  },
);

const warnings = computed(() => app.model.outputs.validationWarnings ?? []);

function clearIfUndef(ref: PlRef | undefined): PlRef | undefined {
  return ref;
}

const binMode = computed(() => app.model.outputs.binMode === true);
</script>

<template>
  <PlSlideModal v-model="open" close-on-outside-click shadow>
    <template #title>Settings</template>

    <PlAlert v-for="(w, i) in warnings" :key="i" :type="w.severity === 'error' ? 'error' : 'warn'">
      {{ w.message }}
    </PlAlert>

    <PlSectionSeparator>Inputs</PlSectionSeparator>
    <PlDropdownRef
      v-model="app.model.args.abundanceRef"
      :options="app.model.outputs.abundanceOptions"
      label="Abundance (reads per clonotype)"
      required
      @update:model-value="
        (v: PlRef | undefined) => (app.model.args.abundanceRef = clearIfUndef(v))
      "
    />
    <PlDropdownRef
      v-model="app.model.args.concentrationColumnRef"
      :options="app.model.outputs.concentrationOptions"
      label="Antigen concentration (per sample)"
      required
    >
      <template #tooltip>
        Per-sample numeric column. The column label is copied into the K_D,app unit annotation —
        prefer a bare unit string (e.g. "nM", "µM") over a phrase.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      v-model="app.model.args.binColumnRef"
      :options="app.model.outputs.binOptions"
      label="FACS bin (optional)"
      clearable
    >
      <template #tooltip>
        Per-sample positive integer. Leave empty to run in no-bin mode — K_D,app values from that
        mode are not comparable to bin-derived results.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      v-model="app.model.args.antigenColumnRef"
      :options="app.model.outputs.antigenOptions"
      label="Antigen (optional)"
      clearable
    />
    <PlTextField
      v-if="app.model.args.antigenColumnRef"
      v-model="app.model.args.targetAntigen"
      label="Target antigen"
      :placeholder="'e.g. Spike-RBD'"
    >
      <template #tooltip>
        Which antigen value to analyse. Required when an antigen column is selected.
      </template>
    </PlTextField>

    <PlAccordionSection label="Read coverage">
      <PlNumberField
        v-model="app.model.args.minReadsPerConcentration"
        label="Min reads per concentration"
        :min-value="1"
        :step="1"
      >
        <template #tooltip>
          Floor applied per clonotype per concentration. Points below the floor are excluded;
          weighting handles remaining coverage variation. Default 3.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.args.minConcentrationPoints"
        label="Min concentration points"
        :min-value="3"
        :step="1"
      >
        <template #tooltip>
          Minimum number of concentration points (after floor filtering) to attempt a fit. Default 5
          — the practical minimum for a 4-parameter Hill fit.
        </template>
      </PlNumberField>
    </PlAccordionSection>

    <PlAccordionSection label="Fit quality">
      <PlNumberField
        v-model="app.model.args.r2ThresholdGood"
        label="R² threshold — Good"
        :min-value="0"
        :max-value="1"
        :step="0.05"
      />
      <PlNumberField
        v-model="app.model.args.r2ThresholdFailed"
        label="R² threshold — Failed"
        :min-value="0"
        :max-value="1"
        :step="0.05"
      />
      <PlNumberField
        v-model="app.model.args.nMin"
        label="Hill coefficient — min"
        :min-value="0"
        :step="0.1"
      />
      <PlNumberField
        v-model="app.model.args.nMax"
        label="Hill coefficient — max"
        :min-value="0"
        :step="0.1"
      >
        <template #tooltip>
          For multimeric antigens (trimeric Spike, dimeric receptors), genuine cooperativity can
          push n into 2–4; raise the max accordingly.
        </template>
      </PlNumberField>
    </PlAccordionSection>

    <PlAccordionSection label="Hook effect">
      <PlNumberField
        v-if="binMode"
        v-model="app.model.args.hookEffectThresholdBin"
        label="Signal drop threshold (bin mode)"
        :min-value="0"
        :step="0.05"
      >
        <template #tooltip>
          Flags non-monotonic signals (potentially genuine tight binders showing a hook effect) when
          the top concentration drops by more than this.
        </template>
      </PlNumberField>
      <PlNumberField
        v-else
        v-model="app.model.args.hookEffectThresholdNoBin"
        label="Signal drop threshold (frequency mode)"
        :min-value="0"
        :step="0.005"
      />
      <PlNumberField
        v-model="app.model.args.hookEffectMinReads"
        label="Min reads for hook check"
        :min-value="0"
        :step="1"
      />
    </PlAccordionSection>

    <PlAccordionSection label="Block label">
      <PlTextField
        v-model="app.model.args.customBlockLabel"
        label="Custom label"
        :placeholder="app.model.args.defaultBlockLabel"
      />
    </PlAccordionSection>
  </PlSlideModal>
</template>
