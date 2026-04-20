<script setup lang="ts">
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

const binMode = computed(() => app.model.outputs.binMode === true);
</script>

<template>
  <PlSlideModal v-model="open" close-on-outside-click shadow>
    <template #title>Inputs &amp; Parameters</template>

    <PlAlert v-for="(w, i) in warnings" :key="i" :type="w.severity === 'error' ? 'error' : 'warn'">
      {{ w.message }}
    </PlAlert>

    <PlSectionSeparator>Inputs</PlSectionSeparator>
    <PlDropdownRef
      v-model="app.model.data.abundanceRef"
      :options="app.model.outputs.abundanceOptions"
      label="Clonotype read counts"
      required
    >
      <template #tooltip>
        Read count per (sample, clonotype) — the MiXCR clonotyping output. Normalized per-sample
        inside the block.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      v-model="app.model.data.concentrationColumnRef"
      :options="app.model.outputs.concentrationOptions"
      label="Antigen concentration"
      required
    >
      <template #tooltip>
        Per-sample numeric column carrying the antigen concentration used for each FACS sort. The
        column label becomes the K_D,app unit — prefer a bare unit (e.g. "nM", "µM") over a phrase.
      </template>
    </PlDropdownRef>
    <PlDropdownRef
      v-model="app.model.data.binColumnRef"
      :options="app.model.outputs.binOptions"
      label="FACS bin"
      clearable
    >
      <template #tooltip>
        Per-sample positive integer identifying the FACS bin. Leave empty to run in no-bin mode —
        K_D,app values from that mode are not comparable to bin-derived results.
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
    <PlTextField
      v-if="app.model.data.antigenColumnRef"
      v-model="app.model.data.targetAntigen"
      label="Target antigen"
      :placeholder="'e.g. Spike-RBD'"
    >
      <template #tooltip>
        Which antigen value to analyse. Required when an antigen column is selected.
      </template>
    </PlTextField>

    <PlAccordionSection label="Read coverage">
      <PlNumberField
        v-model="app.model.data.minReadsPerConcentration"
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
        v-model="app.model.data.minConcentrationPoints"
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
        v-model="app.model.data.r2ThresholdGood"
        label="R² threshold — Good"
        :min-value="0"
        :max-value="1"
        :step="0.05"
      />
      <PlNumberField
        v-model="app.model.data.r2ThresholdFailed"
        label="R² threshold — Failed"
        :min-value="0"
        :max-value="1"
        :step="0.05"
      />
      <PlNumberField
        v-model="app.model.data.nMin"
        label="Hill coefficient — min"
        :min-value="0"
        :step="0.1"
      />
      <PlNumberField
        v-model="app.model.data.nMax"
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
        v-model="app.model.data.hookEffectThresholdBin"
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
        v-model="app.model.data.hookEffectThresholdNoBin"
        label="Signal drop threshold (frequency mode)"
        :min-value="0"
        :step="0.005"
      />
      <PlNumberField
        v-model="app.model.data.hookEffectMinReads"
        label="Min reads for hook check"
        :min-value="0"
        :step="1"
      />
    </PlAccordionSection>

    <PlAccordionSection label="Block label">
      <PlTextField
        v-model="app.model.data.customBlockLabel"
        label="Custom label"
        :placeholder="app.model.data.defaultBlockLabel"
      />
    </PlAccordionSection>
  </PlSlideModal>
</template>
