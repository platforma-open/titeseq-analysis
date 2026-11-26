import type { PColumnSpec } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import { aaCurveDefaults, pinOnVariantAxis, variantCurveDefaults } from "./variantPanel";

const NT_AXIS = {
  name: "pl7.app/variantKey",
  type: "String",
  domain: { "pl7.app/alphabet": "nucleotide" },
} as const;
const REPLICATE_AXIS = {
  name: "pl7.app/vdj/replicate",
  type: "String",
} as const;
const CONCENTRATION_AXIS = {
  name: "pl7.app/vdj/concentration",
  type: "String",
} as const;

const fitted: PColumnSpec = {
  kind: "PColumn",
  name: "pl7.app/vdj/fittedMeanBin",
  valueType: "Double",
  axesSpec: [NT_AXIS, REPLICATE_AXIS, CONCENTRATION_AXIS],
};
const observed: PColumnSpec = { ...fitted, name: "pl7.app/vdj/meanBin" };
const concentration: PColumnSpec = {
  kind: "PColumn",
  name: "pl7.app/vdj/concentrationValue",
  valueType: "Double",
  axesSpec: [CONCENTRATION_AXIS],
};
const aaKey: PColumnSpec = {
  kind: "PColumn",
  name: "pl7.app/vdj/aaVariantKey",
  valueType: "String",
  axesSpec: [NT_AXIS],
};

describe("the variant curve defaults", () => {
  const specs = [observed, fitted, concentration, aaKey];
  const sources = (input: string) =>
    (variantCurveDefaults(specs) ?? [])
      .filter((o) => o.inputName === input)
      .map((o) => o.selectedSource);

  // A scatterplot draws only once every axis of y that x lacks is faceted: grouping and a pinned
  // filter on another column leave the axis free, and the chart asks for a data mapping instead.
  it("facet every axis of the signal that the concentration column lacks", () => {
    const [y] = sources("y");
    const [x] = sources("x");
    const free = (y as PColumnSpec).axesSpec.filter(
      (a) => !(x as PColumnSpec).axesSpec.some((b) => b.name === a.name),
    );
    expect(free.length).toBeGreaterThan(0);
    expect(sources("facetBy")).toEqual(expect.arrayContaining(free));
  });

  it("plot the fitted curve against concentration, the observed points off", () => {
    expect(sources("x")).toEqual([concentration]);
    expect(sources("y")).toEqual([fitted]);
    expect(sources("additionalCurves")).toEqual([fitted]);
  });

  it("colour by replicate", () => {
    expect(sources("grouping")).toEqual([REPLICATE_AXIS]);
  });

  it("facet only the variant and colour nothing on a run without replicates", () => {
    const single = (c: PColumnSpec): PColumnSpec => ({
      ...c,
      axesSpec: [NT_AXIS, CONCENTRATION_AXIS],
    });
    const defaults = variantCurveDefaults([single(observed), single(fitted), concentration]) ?? [];
    expect(defaults.filter((o) => o.inputName === "facetBy").map((o) => o.selectedSource)).toEqual([
      NT_AXIS,
    ]);
    expect(defaults.some((o) => o.inputName === "grouping")).toBe(false);
  });

  it("put every replicate in one panel with the Kd,app marks where the run wrote chart columns", () => {
    const charted = (name: string): PColumnSpec => ({ ...fitted, name });
    const x = charted("pl7.app/vdj/signalConcentration");
    const extent = charted("pl7.app/vdj/curveExtent");
    const mean = charted("pl7.app/vdj/curveCombinedFittedMeanBin");
    const kdLine = charted("pl7.app/vdj/kdLine");
    const kdSe = charted("pl7.app/vdj/kdSeRange");
    expect(variantCurveDefaults([...specs, kdSe, kdLine, mean, extent, x])).toEqual([
      { inputName: "x", selectedSource: x },
      { inputName: "y", selectedSource: extent },
      { inputName: "grouping", selectedSource: REPLICATE_AXIS },
      { inputName: "additionalCurves", selectedSource: fitted },
      { inputName: "additionalCurves", selectedSource: mean },
      { inputName: "additionalCurves", selectedSource: kdLine },
      { inputName: "additionalCurves", selectedSource: kdSe },
    ]);
  });

  it("give no defaults while a column is missing", () => {
    expect(variantCurveDefaults([observed, fitted])).toBeUndefined();
  });
});
describe("the variant axis pin", () => {
  it("matches one variant on the signal's variant axis", () => {
    expect(pinOnVariantAxis([observed, fitted], "pl7.app/vdj/fittedMeanBin", "nt1")).toEqual([
      {
        inputName: "filters",
        selectedSource: NT_AXIS,
        filterType: "equals",
        selectedFilterValues: ["nt1"],
      },
    ]);
  });

  it("pins nothing when the frame has no such signal", () => {
    expect(pinOnVariantAxis([observed], "pl7.app/vdj/fittedMeanBin", "nt1")).toBeUndefined();
  });
});

describe("the amino-acid curves", () => {
  const AA_AXIS = {
    name: "pl7.app/variantKey",
    type: "String",
    domain: { "pl7.app/alphabet": "aminoacid" },
  } as const;
  const onAa = (name: string): PColumnSpec => ({
    kind: "PColumn",
    name,
    valueType: "Double",
    axesSpec: [AA_AXIS, REPLICATE_AXIS, CONCENTRATION_AXIS],
  });
  const aaFitted = onAa("pl7.app/vdj/aaFittedMeanBin");
  const aaMean = onAa("pl7.app/vdj/aaCombinedFittedMeanBin");
  const aaX = onAa("pl7.app/vdj/aaSignalConcentration");

  it("plot the replicates by colour in one panel: each replicate's curve, the mean and the Kd,app marks, on the column that spans them", () => {
    const extent = onAa("pl7.app/vdj/aaCurveExtent");
    const kdLine = onAa("pl7.app/vdj/aaKdLine");
    const kdSe = onAa("pl7.app/vdj/aaKdSeRange");
    const specs = [fitted, aaFitted, aaMean, aaX, extent, kdSe, kdLine];
    expect(aaCurveDefaults(specs)).toEqual([
      { inputName: "x", selectedSource: aaX },
      { inputName: "y", selectedSource: extent },
      { inputName: "grouping", selectedSource: REPLICATE_AXIS },
      { inputName: "additionalCurves", selectedSource: aaFitted },
      { inputName: "additionalCurves", selectedSource: aaMean },
      { inputName: "additionalCurves", selectedSource: kdLine },
      { inputName: "additionalCurves", selectedSource: kdSe },
    ]);
  });

  it("are absent from a run that wrote no amino-acid curves", () => {
    expect(aaCurveDefaults([fitted, observed, concentration])).toBeUndefined();
  });

  it("pin to the amino-acid variant on the amino-acid axis", () => {
    expect(pinOnVariantAxis([fitted, aaFitted], "pl7.app/vdj/aaFittedMeanBin", "aa-1")).toEqual([
      {
        inputName: "filters",
        selectedSource: AA_AXIS,
        filterType: "equals",
        selectedFilterValues: ["aa-1"],
      },
    ]);
  });
});
