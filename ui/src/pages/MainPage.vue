<script setup lang="ts">
import {
  PlAlert,
  PlBlockPage,
  PlDropdownRef,
  PlSectionSeparator,
  PlTextField,
} from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "../app";
import PageHeader from "../components/PageHeader.vue";

const app = useApp();

const warnings = computed(() => app.model.outputs.validationWarnings ?? []);
const hasBlockingError = computed(() => warnings.value.some((w) => w.severity === "error"));
const ready = computed(
  () =>
    app.model.data.abundanceRef !== undefined &&
    app.model.data.concentrationColumnRef !== undefined &&
    !hasBlockingError.value,
);
</script>

<template>
  <PlBlockPage
    v-model:subtitle="app.model.data.customBlockLabel"
    :subtitle-placeholder="app.model.data.defaultBlockLabel"
    title="Overview"
  >
    <template #append>
      <PageHeader />
    </template>

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

    <PlSectionSeparator>Status</PlSectionSeparator>
    <PlAlert v-if="!ready" type="info">
      Select at least an abundance column and a concentration column to run the block. Use
      <strong>Inputs</strong> (top right) to tune fit thresholds and hook-effect parameters.
    </PlAlert>
    <PlAlert v-else type="info">
      Inputs configured. Navigate to a visualization tab (Titration Curves, K_D Distribution,
      Affinity vs Fit Quality, Table) to see results. Use <strong>Inputs</strong> (top right) to
      tune fit thresholds and hook-effect parameters.
    </PlAlert>
  </PlBlockPage>
</template>
