import type { PColumnSpec } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import { variantCurveDefaults, variantCurvePins } from "./variantPanel";

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

describe("the variant curve pins", () => {
  it("match the row key on the curves' own variant axis outside synonymous mode", () => {
    expect(variantCurvePins([fitted, aaKey], "nt-1", false)).toEqual([
      {
        inputName: "filters",
        selectedSource: NT_AXIS,
        filterType: "equals",
        selectedFilterValues: ["nt-1"],
      },
    ]);
  });

  it("match the amino-acid key column on the nucleotide axis in synonymous mode", () => {
    expect(variantCurvePins([fitted, aaKey], "aa-1", true)).toEqual([
      {
        inputName: "filters",
        selectedSource: aaKey,
        filterType: "equals",
        selectedFilterValues: ["aa-1"],
      },
    ]);
  });

  it("pin nothing while no row is open", () => {
    expect(variantCurvePins([fitted, aaKey], undefined, false)).toBeUndefined();
  });

  it("pin nothing when the frame has no fitted curves", () => {
    expect(variantCurvePins([aaKey], "nt-1", false)).toBeUndefined();
  });

  it("pin nothing in synonymous mode without the amino-acid key column", () => {
    expect(variantCurvePins([fitted], "aa-1", true)).toBeUndefined();
  });

  it("ignore an amino-acid key column on another axis", () => {
    const elsewhere: PColumnSpec = {
      ...aaKey,
      axesSpec: [{ name: "pl7.app/other", type: "String" }],
    };
    expect(variantCurvePins([fitted, elsewhere], "aa-1", true)).toBeUndefined();
  });
});

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

  it("plot the observed signal against concentration with the fitted curve", () => {
    expect(sources("x")).toEqual([concentration]);
    expect(sources("y")).toEqual([observed]);
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

  it("give no defaults while a column is missing", () => {
    expect(variantCurveDefaults([observed, fitted])).toBeUndefined();
  });
});
