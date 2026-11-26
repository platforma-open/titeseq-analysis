import type { GraphMakerState, PredefinedGraphOption } from "@milaboratories/graph-maker";
import type { PColumnSpec } from "@platforma-sdk/model";

const OBSERVED_SIGNAL = "pl7.app/vdj/meanBin";
const CONCENTRATION_VALUE = "pl7.app/vdj/concentrationValue";
export const REPLICATE_AXIS = "pl7.app/vdj/replicate";
export const CONCENTRATION_AXIS = "pl7.app/vdj/concentration";

/** The columns of one chart of the variant panel, by role. */
export type CurveChart = {
  /** Each replicate's fitted curve. */
  fitted: string;
  /** With replicates: the mean of the replicate curves. */
  mean: string;
  /** The numeric concentration, on the chart's axes. */
  x: string;
  /** The height of every point of the chart: the column that spans the curves and the marks. */
  extent: string;
  /** A vertical line at Kd,app. */
  kdLine: string;
  /** A bar one standard error either side of Kd,app. */
  kdSe: string;
};

/** The chart of one fitted variant. */
export const VARIANT_CHART: CurveChart = {
  fitted: "pl7.app/vdj/fittedMeanBin",
  mean: "pl7.app/vdj/curveCombinedFittedMeanBin",
  x: "pl7.app/vdj/signalConcentration",
  extent: "pl7.app/vdj/curveExtent",
  kdLine: "pl7.app/vdj/kdLine",
  kdSe: "pl7.app/vdj/kdSeRange",
};

/** Synonymous mode with replicates: the chart of an amino-acid variant, from its nucleotide variants' fits. */
export const AA_CHART: CurveChart = {
  fitted: "pl7.app/vdj/aaFittedMeanBin",
  mean: "pl7.app/vdj/aaCombinedFittedMeanBin",
  x: "pl7.app/vdj/aaSignalConcentration",
  extent: "pl7.app/vdj/aaCurveExtent",
  kdLine: "pl7.app/vdj/aaKdLine",
  kdSe: "pl7.app/vdj/aaKdSeRange",
};

type CurveOption = PredefinedGraphOption<"scatterplot">;

/** The variant panel's curves: points and fitted curves on a log concentration axis. */
export const VARIANT_CURVES_DEFAULT_STATE: GraphMakerState = {
  title: "",
  template: "dots",
  currentTab: null,
  layersSettings: { dots: {} },
  axesSettings: { axisX: { scale: "log" } },
};

/**
 * The variant panel's chart of one fitted variant: observed signal against concentration, with the
 * fitted curve, coloured by replicate. A scatterplot draws only once every axis of y that x lacks
 * is faceted, and the concentration sits on its own axis: so the variant and each replicate get a
 * facet. Undefined while the frame lacks a column the chart needs.
 */
export function variantCurveDefaults(specs: PColumnSpec[]): CurveOption[] | undefined {
  const charted = chartDefaults(specs, VARIANT_CHART);
  if (charted !== undefined) return charted;
  const observed = specs.find((s) => s.name === OBSERVED_SIGNAL);
  const fitted = specs.find((s) => s.name === VARIANT_CHART.fitted);
  const concentration = specs.find((s) => s.name === CONCENTRATION_VALUE);
  if (observed === undefined || fitted === undefined || concentration === undefined) {
    return undefined;
  }
  const free = observed.axesSpec.filter(
    (a) => !concentration.axesSpec.some((c) => c.name === a.name),
  );
  const replicate = observed.axesSpec.find((a) => a.name === REPLICATE_AXIS);
  return [
    { inputName: "x", selectedSource: concentration },
    { inputName: "y", selectedSource: fitted },
    ...free.map((a): CurveOption => ({ inputName: "facetBy", selectedSource: a })),
    ...(replicate ? [{ inputName: "grouping" as const, selectedSource: replicate }] : []),
    { inputName: "additionalCurves", selectedSource: fitted },
  ];
}

/**
 * A chart where the run wrote its columns: every replicate's curve in one panel, coloured by
 * replicate, the mean curve over replicates, and Kd,app with its standard error. Every column
 * shares one set of axes, so nothing is faceted. A variant without Kd,app has no marks, but the
 * columns are always written. Undefined on a run made before the columns.
 */
function chartDefaults(specs: PColumnSpec[], chart: CurveChart): CurveOption[] | undefined {
  const named = (name: string) => specs.find((s) => s.name === name);
  const fitted = named(chart.fitted);
  const extent = named(chart.extent);
  const concentration = named(chart.x);
  if (!fitted || !extent || !concentration) return undefined;
  const replicate = fitted.axesSpec.find((a) => a.name === REPLICATE_AXIS);
  const curves = [fitted, named(chart.mean), named(chart.kdLine), named(chart.kdSe)];
  return [
    { inputName: "x", selectedSource: concentration },
    // The points are hidden: y is the column that spans the curves.
    { inputName: "y", selectedSource: extent },
    ...(replicate ? [{ inputName: "grouping" as const, selectedSource: replicate }] : []),
    ...curves
      .filter((s) => s !== undefined)
      .map((s): CurveOption => ({ inputName: "additionalCurves", selectedSource: s })),
  ];
}

/**
 * Pins a chart to one variant on the variant axis of the column `signalName`: the axis that is
 * neither the replicate nor the concentration. A hidden graph-maker filter: applied to the chart's
 * data, re-applied when it changes, never saved with the chart. Undefined without the column.
 */
export function pinOnVariantAxis(
  specs: PColumnSpec[],
  signalName: string,
  key: string,
): CurveOption[] | undefined {
  const axis = specs
    .find((s) => s.name === signalName)
    ?.axesSpec.find((a) => a.name !== REPLICATE_AXIS && a.name !== CONCENTRATION_AXIS);
  if (axis === undefined) return undefined;
  return [
    {
      inputName: "filters",
      selectedSource: axis,
      filterType: "equals",
      selectedFilterValues: [key],
    },
  ];
}

/**
 * Synonymous mode with replicates: the panel's chart of an amino-acid variant, made from its
 * nucleotide variants' fits. Undefined on a run that wrote no amino-acid curves.
 */
export function aaCurveDefaults(specs: PColumnSpec[]): CurveOption[] | undefined {
  return chartDefaults(specs, AA_CHART);
}
