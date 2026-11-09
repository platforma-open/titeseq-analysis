import type { GraphMakerState } from "@milaboratories/graph-maker";
// The package root exports only `GraphMaker` and types; this runtime function is reachable only
// through the `./*` wildcard, which maps literally and so needs the `.js`.
import { createDefaultMapping } from "@milaboratories/graph-maker/dist/dataBindAes.js";
import type { AxisSpec } from "@platforma-sdk/model";
import type { Ref } from "vue";
import { watch } from "vue";
import { useApp } from "../app";

/**
 * The chart's key for an axis, byte for byte: keys sorted, `domain` always present (`{}` when the
 * spec has none).
 */
export function axisId(axis: AxisSpec): string {
  const domain = Object.fromEntries(
    Object.entries(axis.domain ?? {}).sort(([a], [b]) => a.localeCompare(b)),
  );
  return JSON.stringify({ domain, kind: "axis", name: axis.name, type: axis.type });
}

/**
 * Keeps a category chart's concentrations in numeric order.
 *
 * The String concentration axis sorts as text (1e-6 before 3.16e-6 before 1e-7), and the chart
 * refuses the numeric concentrationValue column as a grouping, so the order lives in
 * `dataBindAes`. The chart owns that entry: it deletes it while the frame has no values (a re-run)
 * and rebuilds it in text order when they return. So the order is corrected whenever the chart's
 * entry differs, never created ahead of it. The entry must be complete: the chart reads
 * `mapping[value].aes` for every value in `order`.
 *
 * @param state the chart's stored state; written only to correct the order
 * @param axis the concentration axis the chart groups by
 */
export function useConcentrationOrder(
  state: Ref<GraphMakerState | undefined>,
  axis: Ref<AxisSpec | undefined>,
) {
  const app = useApp();
  watch(
    [
      axis,
      () => app.model.outputs.concentrationOrder,
      () => axis.value && state.value?.dataBindAes?.[axisId(axis.value)],
    ],
    ([concentration, numeric, entry]) => {
      const current = state.value;
      if (!concentration || !numeric?.length || !current || entry?.type !== "categorical") return;
      const present = new Set(entry.order);
      const order = [
        ...numeric.filter((v) => present.has(v)),
        ...entry.order.filter((v) => !numeric.includes(v)),
      ];
      if (order.every((v, i) => v === entry.order[i])) return;
      state.value = {
        ...current,
        dataBindAes: {
          ...current.dataBindAes,
          [axisId(concentration)]: createDefaultMapping(order, "bright"),
        },
      };
    },
    { immediate: true },
  );
}
