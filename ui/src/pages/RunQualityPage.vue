<script setup lang="ts">
import type { PlChartStackedBarSettings } from "@platforma-sdk/ui-vue";
import { PlChartStackedBar, PlTabs } from "@platforma-sdk/ui-vue";
import { computed, ref } from "vue";
import { useApp } from "../app";
import CountsGraph from "../components/CountsGraph.vue";
import RunQualityGate from "../components/RunQualityGate.vue";
import TiteseqPage from "../components/TiteseqPage.vue";

const app = useApp();

type Tab = "outcomes" | "failed" | "raw";
const tabs: { label: string; value: Tab }[] = [
  { label: "Fit Outcomes", value: "outcomes" },
  { label: "Failed Fits", value: "failed" },
  { label: "Raw Counts", value: "raw" },
];
const tab = ref<Tab>("outcomes");

// Best to worst, so every bar reads left to right the same way. The colors are the ones the
// workflow declares for the GraphMaker charts: its "triadic" palette.
const CLASSES = [
  { name: "Good", color: "#99E099", description: "Kd,app measured, with a small error." },
  { name: "Partial", color: "#FFCB8F", description: "Kd,app measured, with a wide error." },
  { name: "Above range", color: "#99CCFF", description: "Weaker than the titration can measure." },
  { name: "Below range", color: "#C1ADFF", description: "Tighter than the titration can measure." },
  { name: "Failed", color: "#F05670", description: "No Kd,app." },
];
// Every reason a single fit can report, as the workflow declares them.
const REASONS = [
  { name: "Insufficient points", color: "#F05670" },
  { name: "Insufficient reads", color: "#FFCB8F" },
  { name: "Bright without antigen", color: "#C1ADFF" },
  { name: "Non monotonic signal", color: "#99CCFF" },
  { name: "Incomplete titration", color: "#90E0E0" },
  { name: "Convergence failure", color: "#FAAAFA" },
  { name: "Low R²", color: "#CBEB67" },
  { name: "High Kd,app SE", color: "#FFADBA" },
  { name: "Hill coefficient out of range", color: "#D3D7E0" },
];
const OTHER_REASON_COLOR = "#929BAD";

type Segment = PlChartStackedBarSettings["data"][number];
type Chart = { key: string; empty: boolean; settings: PlChartStackedBarSettings };

const replicateName = (replicate: string) =>
  replicate === app.model.outputs.runQuality?.allReplicates ? replicate : `Replicate ${replicate}`;

const total = (segments: Segment[]) => segments.reduce((sum, s) => sum + s.value, 0);

// Every replicate of the run, from the raw counts: a replicate with no counted fit has no other row.
const replicates = computed(() => [
  ...new Set((app.model.outputs.runQuality?.raw ?? []).map((r) => r.replicate)),
]);

/**
 * All replicates together, then one chart per replicate; one chart for a run without replicates.
 * Only the first chart carries the legend: it counts every replicate, so it has every segment.
 */
function charts<R extends { replicate?: string }>(
  rows: R[],
  title: (count: number) => string,
  segments: (rows: R[]) => Segment[],
): Chart[] {
  return replicates.value.map((replicate, i) => {
    const data = segments(rows.filter((r) => r.replicate === replicate));
    const prefix = replicate === undefined ? "" : `${replicateName(replicate)}: `;
    return {
      key: replicate ?? "",
      empty: data.length === 0,
      settings: {
        title: prefix + title(total(data)),
        data,
        maxLegendsInColumn: 5,
        showLegends: i === 0,
      },
    };
  });
}

const outcomeCharts = computed(() =>
  charts(
    app.model.outputs.runQuality?.outcomes ?? [],
    (n) => `${n.toLocaleString()} fits`,
    (rows) =>
      CLASSES.flatMap((c) => {
        const row = rows.find((r) => r.affinityClass === c.name);
        return row
          ? [{ label: c.name, value: row.count, color: c.color, description: c.description }]
          : [];
      }),
  ),
);

// The declared reasons in their order, then any other the run reports.
const reasonColors = computed(() => {
  const colors = new Map(REASONS.map((r) => [r.name, r.color]));
  for (const r of app.model.outputs.runQuality?.failureReasons ?? [])
    if (!colors.has(r.failureReason)) colors.set(r.failureReason, OTHER_REASON_COLOR);
  return colors;
});

const reasonCharts = computed(() =>
  charts(
    app.model.outputs.runQuality?.failureReasons ?? [],
    (n) => `${n.toLocaleString()} failed fits`,
    (rows) =>
      [...reasonColors.value.entries()].flatMap(([reason, color]) => {
        const row = rows.find((r) => r.failureReason === reason);
        return row ? [{ label: reason, value: row.count, color }] : [];
      }),
  ),
);

/**
 * The two tabs of counts: the bars, then the GraphMaker chart of `column`. `none` is what a chart
 * with nothing to draw says; `blank` leaves the tab with that sentence alone.
 */
const countTabs = computed(() => [
  {
    tab: "outcomes" as const,
    charts: outcomeCharts.value,
    none: "No fit is counted.",
    blank: false,
    state: "graphStateFitOutcomes" as const,
    column: "pl7.app/vdj/fitOutcomeCount",
  },
  {
    tab: "failed" as const,
    charts: reasonCharts.value,
    none: "No fit failed.",
    blank: (app.model.outputs.runQuality?.failureReasons ?? []).length === 0,
    state: "graphStateFailureReasons" as const,
    column: "pl7.app/vdj/fitFailureCount",
  },
]);

/** Every fit: a row per class and failure reason, a column per replicate. */
const rawTable = computed(() => {
  const quality = app.model.outputs.runQuality;
  const raw = quality?.raw ?? [];
  const leftOut = quality?.leftOutReasons ?? [];
  const order = (r: { affinityClass: string }) =>
    CLASSES.findIndex((c) => c.name === r.affinityClass);
  const label = (r: { affinityClass: string; failureReason: string | null }) =>
    r.failureReason === null ? r.affinityClass : `${r.affinityClass}: ${r.failureReason}`;
  const labels = [...new Set([...raw].sort((a, b) => order(a) - order(b)).map(label))];
  const count = (name: string, replicate: string | undefined) =>
    raw
      .filter((r) => label(r) === name && r.replicate === replicate)
      .reduce((sum, r) => sum + r.count, 0);
  const rows = labels.map((name) => ({
    name,
    excluded: leftOut.some((reason) => name === `Failed: ${reason}`),
    counts: replicates.value.map((rep) => count(name, rep)),
  }));
  const totals = replicates.value.map((_, i) => rows.reduce((sum, row) => sum + row.counts[i], 0));
  return {
    headers: replicates.value.map((rep) => (rep === undefined ? "Fits" : replicateName(rep))),
    note: leftOut.length > 0 ? "The graphs leave out the rows marked *." : undefined,
    rows,
    totals,
  };
});
</script>

<template>
  <TiteseqPage title="Run Quality" mode="table">
    <RunQualityGate>
      <div class="run-quality">
        <PlTabs v-model="tab" :options="tabs" />
        <template v-for="t in countTabs" :key="t.tab">
          <section v-if="tab === t.tab">
            <p v-if="t.blank">{{ t.none }}</p>
            <template v-else>
              <template v-for="c in t.charts" :key="c.key">
                <div v-if="c.empty" class="run-quality-empty">
                  <div class="run-quality-empty__title">{{ c.settings.title }}</div>
                  <p>{{ t.none }}</p>
                </div>
                <PlChartStackedBar v-else :settings="c.settings" />
              </template>
              <div class="run-quality-graph">
                <CountsGraph v-model="app.model.data[t.state]" :column-name="t.column" />
              </div>
            </template>
          </section>
        </template>
        <section v-if="tab === 'raw'">
          <p v-if="rawTable.note">{{ rawTable.note }}</p>
          <table class="raw-counts">
            <thead>
              <tr>
                <th>Outcome</th>
                <th v-for="h in rawTable.headers" :key="h">{{ h }}</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="row in rawTable.rows" :key="row.name">
                <td>{{ row.name }}{{ row.excluded ? " *" : "" }}</td>
                <td v-for="(n, i) in row.counts" :key="i">{{ n.toLocaleString() }}</td>
              </tr>
            </tbody>
            <tfoot>
              <tr>
                <td>Total</td>
                <td v-for="(n, i) in rawTable.totals" :key="i">{{ n.toLocaleString() }}</td>
              </tr>
            </tfoot>
          </table>
        </section>
      </div>
    </RunQualityGate>
  </TiteseqPage>
</template>

<style scoped>
.run-quality {
  display: flex;
  flex-direction: column;
  gap: 32px;
  max-width: 960px;
}
.run-quality section {
  display: flex;
  flex-direction: column;
  gap: 16px;
}
.run-quality-graph {
  height: 640px;
}
/* The title of PlChartStackedBar, for a replicate with nothing to draw. */
.run-quality-empty__title {
  font-size: 20px;
  font-weight: 500;
  line-height: 24px;
  letter-spacing: -0.2px;
  margin-bottom: 8px;
}
.run-quality p {
  margin: 0;
}
.raw-counts {
  border-collapse: collapse;
  align-self: flex-start;
}
.raw-counts th,
.raw-counts td {
  padding: 6px 16px 6px 0;
  text-align: left;
  border-bottom: 1px solid var(--border-color-div-grey);
}
.raw-counts th:not(:first-child),
.raw-counts td:not(:first-child) {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
.raw-counts tfoot td {
  font-weight: 600;
}
</style>
