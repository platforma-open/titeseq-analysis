import type { PColumnSpec } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import { variantCurvePins } from "./variantPanel";

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
