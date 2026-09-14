import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import { kind } from "./index";

const parse = (v: unknown) => kind.parseInitializationParams(v);

const ref = (name: string): PlRef => ({ __isRef: true, blockId: "b1", name });

const COLUMN_REFS = [
  "abundanceRef",
  "concentrationColumnRef",
  "binColumnRef",
  "antigenColumnRef",
  "sortFractionColumnRef",
];

describe("the envelope", () => {
  it("accepts an empty object — a block may be created with nothing pinned", () => {
    expect(parse({})).toEqual({});
  });

  it.each([undefined, null, 42, "params", [], true])("refuses %o as a params object", (v) => {
    expect(() => parse(v)).toThrow();
  });

  it("drops keys the contract does not name", () => {
    expect(parse({ nMin: 1, defaultBlockLabel: "x", settingsOpen: true, tableState: {} })).toEqual({
      nMin: 1,
    });
  });
});

describe("the column references", () => {
  it.each(COLUMN_REFS)("accepts a reference for %s", (field) => {
    expect(parse({ [field]: ref("col") })).toEqual({ [field]: ref("col") });
  });

  it.each(COLUMN_REFS)("refuses a bare object for %s", (field) => {
    expect(() => parse({ [field]: { blockId: "b1", name: "col" } })).toThrow(`'${field}' must be`);
  });

  it.each([null, "b1/col", 7, []])("refuses %o as a reference", (v) => {
    expect(() => parse({ abundanceRef: v })).toThrow("'abundanceRef' must be");
  });
});

describe("the bounded settings", () => {
  it.each([1, 3, 100])("accepts %o reads per concentration", (v) => {
    expect(parse({ minReadsPerConcentration: v })).toEqual({ minReadsPerConcentration: v });
  });

  it.each([0, 0.9, -1])("refuses %o reads per concentration", (v) => {
    expect(() => parse({ minReadsPerConcentration: v })).toThrow(
      "'minReadsPerConcentration' must be",
    );
  });

  it.each([3, 5, 12])("accepts %o concentration points", (v) => {
    expect(parse({ minConcentrationPoints: v })).toEqual({ minConcentrationPoints: v });
  });

  it.each([0, 2, 2.99])("refuses %o concentration points", (v) => {
    expect(() => parse({ minConcentrationPoints: v })).toThrow("'minConcentrationPoints' must be");
  });

  it.each(["r2ThresholdGood", "r2ThresholdFailed"])("accepts 0, 0.5 and 1 for %s", (field) => {
    for (const v of [0, 0.5, 1]) expect(parse({ [field]: v })).toEqual({ [field]: v });
  });

  it.each(["r2ThresholdGood", "r2ThresholdFailed"])("refuses a value outside 0..1 for %s", (f) => {
    for (const v of [-0.01, 1.01, 2]) expect(() => parse({ [f]: v })).toThrow(`'${f}' must be`);
  });
});

describe("the unbounded settings", () => {
  const FREE = [
    "nMin",
    "nMax",
    "hookEffectThresholdBin",
    "hookEffectThresholdNoBin",
    "hookEffectMinReads",
  ];

  it.each(FREE)("accepts any real number for %s", (field) => {
    for (const v of [0, 0.02, 2, 1000]) expect(parse({ [field]: v })).toEqual({ [field]: v });
  });

  it.each(FREE)("refuses NaN and the infinities for %s", (field) => {
    for (const v of [Number.NaN, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY]) {
      expect(() => parse({ [field]: v })).toThrow(`'${field}' must be`);
    }
  });

  it.each(["2", null, {}, []])("refuses %o as a number", (v) => {
    expect(() => parse({ nMin: v })).toThrow("'nMin' must be");
  });
});

describe("the cross-field rules the block enforces at run time", () => {
  // Each of these is a state the settings drawer leaves behind mid-edit, and `.args()` refuses
  // it when the block runs. Refusing it here would mean a block reachable by hand cannot be
  // carried by a template.
  it("accepts a Hill range the wrong way round", () => {
    expect(parse({ nMin: 2, nMax: 0.5 })).toEqual({ nMin: 2, nMax: 0.5 });
  });

  it("accepts a failed R² threshold above the good one", () => {
    expect(parse({ r2ThresholdGood: 0.5, r2ThresholdFailed: 0.8 })).toEqual({
      r2ThresholdGood: 0.5,
      r2ThresholdFailed: 0.8,
    });
  });

  it("accepts an antigen column with no target antigen named", () => {
    expect(parse({ antigenColumnRef: ref("antigen") })).toEqual({
      antigenColumnRef: ref("antigen"),
    });
  });

  it("accepts a sort fraction column with no bin column", () => {
    expect(parse({ sortFractionColumnRef: ref("sf") })).toEqual({
      sortFractionColumnRef: ref("sf"),
    });
  });
});

describe("the text settings", () => {
  it("accepts a target antigen and a custom label, empty strings included", () => {
    expect(parse({ targetAntigen: "", customBlockLabel: "" })).toEqual({
      targetAntigen: "",
      customBlockLabel: "",
    });
    expect(parse({ targetAntigen: "CD19", customBlockLabel: "run 4" })).toEqual({
      targetAntigen: "CD19",
      customBlockLabel: "run 4",
    });
  });

  it.each([7, null, {}])("refuses %o as a target antigen", (v) => {
    expect(() => parse({ targetAntigen: v })).toThrow("'targetAntigen' must be");
  });
});

describe("identity", () => {
  it("names this package and its published version", () => {
    expect(kind.name).toBe("@platforma-open/platforma-open.titeseq-analysis.kind");
    expect(kind.version).toMatch(/^\d+\.\d+\.\d+/);
  });
});
