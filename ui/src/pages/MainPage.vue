<script setup lang="ts">
import { PlAlert, PlBlockPage, PlTextField } from "@platforma-sdk/ui-vue";
import { useApp } from "../app";

const app = useApp();
</script>

<template>
  <PlBlockPage>
    <PlTextField v-model="app.model.data.name" label="Enter your name" :clearable="() => ''" />

    <PlAlert v-if="app.model.outputs.tengoMessage" type="success">
      {{ app.model.outputs.tengoMessage }}
    </PlAlert>

    <PlAlert v-if="app.model.outputs.pythonMessage" type="success">
      {{ app.model.outputs.pythonMessage }}
    </PlAlert>
  </PlBlockPage>
</template>
