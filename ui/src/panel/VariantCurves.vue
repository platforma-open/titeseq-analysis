<script setup lang="ts">
/**
 * The variant panel's Curves tab: the opened row's titration curves.
 *
 * - One fitted variant: its observed points and fitted curve, a facet per replicate.
 * - Synonymous mode: the row is an amino-acid variant and the fits are its synonymous variants', so
 *   each synonymous variant has a tab.
 * - Synonymous mode with replicates: the amino-acid variant's own chart. Each replicate's points
 *   and curve, averaged over its synonymous variants, and the mean curve over replicates.
 */
import {
  AA_CHART,
  REPLICATE_AXIS,
  VARIANT_CHART,
  VARIANT_CURVES_DEFAULT_STATE,
  aaCurveDefaults,
  pinOnVariantAxis,
  variantCurveDefaults,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import { AxisId, ColumnId } from "@milaboratories/pf-plots";
import type { PColumnIdAndSpec } from "@platforma-sdk/model";
import { PlTabs } from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
import { useApp } from "../app";
import type { CurveChart } from "@platforma-open/platforma-open.titeseq-analysis.model";
import type { PanelVariant } from "./panelVariants";
import { AA_VARIANT_KEY, handleOf, hasCurve, loadSynonymousVariants } from "./panelVariants";

const props = defineProps<{ rowKey: string }>();
// `true` once the open row is known to have no curve to draw.
const emit = defineEmits<{ empty: [boolean] }>();

const app = useApp();

const INK = "#110529";

const synonymous = computed(() => app.model.outputs.runShape?.synonymous === true);
const replicates = computed(() => app.model.outputs.runShape?.replicates === true);
const cols = computed(() => app.model.outputs.titrationCurvesPfCols);
const specs = computed(() => cols.value?.map((c) => c.spec) ?? []);

// The amino-acid chart, where the run wrote amino-acid curves; a run made before them falls back
// to the synonymous-variant tabs.
const aaDefaults = computed(() =>
  synonymous.value && replicates.value ? aaCurveDefaults(specs.value) : undefined,
);
const aaMode = computed(() => aaDefaults.value !== undefined);

// A synonymous-mode run made before the block linked each fitted variant to its amino-acid variant:
// its curves cannot be found for a row. A run with the amino-acid chart needs no link.
const unlinkedRun = computed(
  () =>
    synonymous.value &&
    !aaMode.value &&
    cols.value !== undefined &&
    !cols.value.some((c) => c.spec.name === AA_VARIANT_KEY),
);

// The fitted variants behind the open row. `undefined` while loading.
const variants = ref<PanelVariant[] | undefined>(undefined);
const picked = ref<string | undefined>(undefined);
const loadError = ref<string | undefined>(undefined);

// Fires only when a source's value changes, not on every output patch.
watch(
  [
    () => props.rowKey,
    () => handleOf(app.model.outputs.titrationCurvesPf),
    synonymous,
    aaMode,
    // The curve check reads the frame's columns, which can arrive after its handle.
    () => cols.value !== undefined,
  ],
  async ([key, frame, synonymous, aa, listed], _old, onCleanup) => {
    // A later firing supersedes this load.
    let stale = false;
    onCleanup(() => {
      stale = true;
    });
    variants.value = undefined;
    picked.value = undefined;
    loadError.value = undefined;
    if (!frame || !listed) return;
    try {
      // In synonymous mode without the amino-acid chart the fits are the row's synonymous
      // variants'. Otherwise the row is charted itself, where it has a curve.
      let found: PanelVariant[];
      if (synonymous && !aa) {
        found = await loadSynonymousVariants(frame, key);
      } else {
        const signal = (aa ? AA_CHART : VARIANT_CHART).fitted;
        found = (await hasCurve(frame, cols.value ?? [], signal, key)) ? [{ key, label: key }] : [];
      }
      if (stale) return;
      variants.value = found;
      picked.value = found[0]?.key;
      emit("empty", found.length === 0);
    } catch (error) {
      if (!stale) loadError.value = error instanceof Error ? error.message : String(error);
    }
  },
  { immediate: true },
);

const pickedVariant = computed(() => variants.value?.find((v) => v.key === picked.value));
const variantTabs = computed(() =>
  (variants.value ?? []).map((v) => ({ value: v.key, label: v.label })),
);

const heading = computed(() => {
  if (aaMode.value) return "Replicate curves of the amino acid variant";
  const base = replicates.value ? "Replicate curves" : "Titration curve";
  return synonymous.value && pickedVariant.value ? `${base} of ${pickedVariant.value.label}` : base;
});

// A chart input as graph-maker keys it in its state.
const curveKey = (col: PColumnIdAndSpec) =>
  new ColumnId({ name: col.columnId, type: col.spec.valueType }).toCanonicalString();

const SOLID = { lineShape: "solid" as const, opacity: 1, showDots: false };

/** The colour of the column's replicate; without the replicate axis, ink. */
function replicateColor(col: PColumnIdAndSpec) {
  const replicate = col.spec.axesSpec.find((a) => a.name === REPLICATE_AXIS);
  return replicate
    ? { type: "grouping" as const, value: AxisId.fromAxisSpec(replicate).toCanonicalString() }
    : INK;
}

/**
 * The chart's state as the block styles it: curves carry the information, so the points are hidden
 * and there is no grid. With replicates each fitted curve is a dimmed dashed line in its
 * replicate's colour, under the mean over replicates as a solid line in ink; without them the
 * fitted curve is that solid line. Kd,app is a dot on the mean curve, each replicate's own Kd,app
 * a dot on its curve in its colour, and the standard error a faint band: one line the run draws
 * back and forth across it, wide enough here to close the gaps.
 */
function styledCurvesState(chart: CurveChart, cols: PColumnIdAndSpec[]): GraphMakerState {
  const named = (name: string) => cols.find((c) => c.spec.name === name);
  const fitted = named(chart.fitted);
  const mean = named(chart.mean);
  const kdPoint = named(chart.kdPoint);
  const kdReplicatePoint = named(chart.kdReplicatePoint);
  const kdSe = named(chart.kdSe);
  const replicated = fitted?.spec.axesSpec.some((a) => a.name === REPLICATE_AXIS) === true;
  const perCurve = {
    ...(fitted && {
      [curveKey(fitted)]: replicated
        ? {
            ...SOLID,
            lineShape: "dashed" as const,
            lineWidth: 1.5,
            lineColor: replicateColor(fitted),
            opacity: 0.7,
          }
        : { ...SOLID, lineWidth: 2, lineColor: INK },
    }),
    ...(mean && { [curveKey(mean)]: { ...SOLID, lineWidth: 2.5, lineColor: INK } }),
    ...(kdPoint && {
      [curveKey(kdPoint)]: { ...SOLID, lineWidth: 1, lineColor: INK, showDots: true },
    }),
    ...(kdReplicatePoint && {
      [curveKey(kdReplicatePoint)]: {
        ...SOLID,
        lineWidth: 1,
        lineColor: replicateColor(kdReplicatePoint),
        showDots: true,
      },
    }),
    ...(kdSe && { [curveKey(kdSe)]: { ...SOLID, lineWidth: 6, lineColor: INK, opacity: 0.15 } }),
  };
  return {
    ...VARIANT_CURVES_DEFAULT_STATE,
    layersSettings: {
      dots: { dotSize: 0 },
      ...(Object.keys(perCurve).length > 0 && { additionalCurves: { perCurve } }),
    },
    axesSettings: {
      axisX: { scale: "log", gridlines: false },
      axisY: { gridlines: false },
      other: { frame: "left-bottom" },
    },
  };
}

const chart = computed(() => (aaMode.value ? AA_CHART : VARIANT_CHART));

// One saved state per run shape: a run with replicates colours its curves by replicate, one
// without has no replicate to colour by, so neither state fits the other.
const stateField = computed(() => {
  if (aaMode.value) return "graphStateVariantAaCurves";
  return replicates.value ? "graphStateVariantReplicateCurves" : "graphStateVariantCurves";
});

// The curve styles are keyed by column id, so a shape's state is styled when the frame's columns
// first name its curves. A re-run can add curves (a run made before the chart columns has the
// fitted curve only), so a saved state gets the style of each curve it has none for. A style the
// state already holds is the user's and stays. The chart reads its state only when it mounts, so
// `styledCount` remounts it.
const styledCount = ref(0);
watch(
  [cols, stateField, chart],
  ([cols, field, chart]) => {
    if (cols === undefined) return;
    const styled = styledCurvesState(chart, cols);
    const saved = app.model.data[field];
    if (saved === undefined) {
      app.model.data[field] = styled;
      return;
    }
    const kept = saved.layersSettings?.additionalCurves?.perCurve ?? {};
    const missing = Object.entries(styled.layersSettings?.additionalCurves?.perCurve ?? {}).filter(
      ([key]) => kept[key] === undefined,
    );
    if (missing.length === 0) return;
    app.model.data[field] = {
      ...saved,
      layersSettings: {
        ...saved.layersSettings,
        additionalCurves: {
          ...saved.layersSettings?.additionalCurves,
          perCurve: { ...kept, ...Object.fromEntries(missing) },
        },
      },
    };
    styledCount.value++;
  },
  { immediate: true },
);

const curvesState = computed({
  get: (): GraphMakerState => app.model.data[stateField.value] ?? VARIANT_CURVES_DEFAULT_STATE,
  set: (state: GraphMakerState) => {
    app.model.data[stateField.value] = state;
  },
});

const curveDefaults = computed((): PredefinedGraphOption<"scatterplot">[] | undefined =>
  aaMode.value ? aaDefaults.value : variantCurveDefaults(specs.value),
);

const curvePins = computed(() => {
  if (aaMode.value) return pinOnVariantAxis(specs.value, AA_CHART.fitted, props.rowKey);
  return picked.value === undefined
    ? undefined
    : pinOnVariantAxis(specs.value, VARIANT_CHART.fitted, picked.value);
});

const chartKey = computed(
  () => `${stateField.value}-${aaMode.value ? props.rowKey : picked.value}-${styledCount.value}`,
);
// The state must exist in data before the chart binds it, or the chart never settles.
const ready = computed(
  () =>
    (variants.value?.length ?? 0) > 0 &&
    curvePins.value !== undefined &&
    app.model.data[stateField.value] !== undefined,
);
</script>

<template>
  <div class="titeseq-variant-curves">
    <div v-if="loadError" class="titeseq-variant-curves-note">
      The variant's curves did not load: {{ loadError }}
    </div>
    <div v-else-if="variants === undefined" class="titeseq-variant-curves-note">Loading…</div>
    <div v-else-if="unlinkedRun" class="titeseq-variant-curves-note">
      These results were made before the block showed a variant's curves here. Run the block again
      to see them.
    </div>
    <div v-else-if="!ready" class="titeseq-variant-curves-note">
      No titration curves for this variant.
    </div>
    <template v-else>
      <PlTabs v-if="!aaMode && variantTabs.length > 1" v-model="picked" :options="variantTabs" />
      <div class="titeseq-variant-curves-heading">{{ heading }}</div>
      <div class="titeseq-variant-curves-chart">
        <GraphMaker
          :key="chartKey"
          v-model="curvesState"
          chart-type="scatterplot"
          :p-frame="app.model.outputs.titrationCurvesPf"
          :default-options="curveDefaults"
          :fixed-options="curvePins"
          default-palette="bright"
        />
      </div>
    </template>
  </div>
</template>

<style scoped>
.titeseq-variant-curves {
  flex: 1;
  display: flex;
  flex-direction: column;
  gap: 12px;
  min-height: 480px;
}
.titeseq-variant-curves-heading {
  font-size: 13px;
  font-weight: 600;
}
/* A chart fills its container, and a container with no height collapses to nothing. */
.titeseq-variant-curves-chart {
  height: 620px;
  display: flex;
  flex-direction: column;
}
.titeseq-variant-curves-note {
  padding: 16px;
  color: var(--txt-03);
}
</style>
