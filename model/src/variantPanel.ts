import type { AxisSpec, PColumnSpec } from "@platforma-sdk/model";

const OBSERVED_SIGNAL = "pl7.app/vdj/meanBin";
const FITTED_SIGNAL = "pl7.app/vdj/fittedMeanBin";
const CONCENTRATION_VALUE = "pl7.app/vdj/concentrationValue";
const AA_VARIANT_KEY = "pl7.app/vdj/aaVariantKey";
const REPLICATE_AXIS = "pl7.app/vdj/replicate";
const CONCENTRATION_AXIS = "pl7.app/vdj/concentration";

/**
 * A hidden graph-maker filter, in the shape of its `fixedOptions` entries. Graph-maker applies it
 * to the chart's data, re-applies it when it changes, and never saves it with the chart.
 */
export type ChartPin = {
  inputName: "filters";
  selectedSource: AxisSpec | PColumnSpec;
  filterType: "equals";
  selectedFilterValues: string[];
};

/** A graph-maker default option, in the shape of its `defaultOptions` entries. */
export type ChartDefault = {
  inputName: "x" | "y" | "facetBy" | "grouping" | "additionalCurves";
  selectedSource: AxisSpec | PColumnSpec;
};

/**
 * The variant panel's chart: observed signal against concentration, with the fitted curve. A
 * scatterplot draws only once every axis of y that x lacks is faceted; grouping and a pinned
 * filter on another column leave an axis free. So each variant and each replicate gets a facet,
 * and colour marks the replicate. Undefined while the frame lacks a column the chart needs.
 */
export function variantCurveDefaults(specs: PColumnSpec[]): ChartDefault[] | undefined {
  const observed = specs.find((s) => s.name === OBSERVED_SIGNAL);
  const fitted = specs.find((s) => s.name === FITTED_SIGNAL);
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
    { inputName: "y", selectedSource: observed },
    ...free.map((a): ChartDefault => ({ inputName: "facetBy", selectedSource: a })),
    ...(replicate ? [{ inputName: "grouping" as const, selectedSource: replicate }] : []),
    { inputName: "additionalCurves", selectedSource: fitted },
  ];
}

/**
 * Pins the titration curves to the variant a table row names. Outside synonymous mode the row key
 * is the curves' own variant key. In synonymous mode the row is an amino-acid variant and the
 * curves are keyed by nucleotide variant, so the pin matches the amino-acid key column on the
 * nucleotide axis: every synonymous variant of that row. Undefined when nothing is open or the
 * frame lacks the column the pin needs.
 */
export function variantCurvePins(
  specs: PColumnSpec[],
  key: string | undefined,
  synonymous: boolean,
): ChartPin[] | undefined {
  if (key === undefined) return undefined;
  const signal = specs.find((s) => s.name === FITTED_SIGNAL);
  if (signal === undefined) return undefined;
  const variantAxis = signal.axesSpec.find(
    (a) => a.name !== REPLICATE_AXIS && a.name !== CONCENTRATION_AXIS,
  );
  if (variantAxis === undefined) return undefined;
  const source = synonymous
    ? specs.find(
        (s) =>
          s.name === AA_VARIANT_KEY &&
          s.axesSpec.length === 1 &&
          s.axesSpec[0].name === variantAxis.name,
      )
    : variantAxis;
  if (source === undefined) return undefined;
  return [
    {
      inputName: "filters",
      selectedSource: source,
      filterType: "equals",
      selectedFilterValues: [key],
    },
  ];
}
