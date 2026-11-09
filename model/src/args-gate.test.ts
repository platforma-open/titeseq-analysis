import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockData } from "./index";
import { deriveArgs, initBlockData } from "./index";

const ref = (name: string): PlRef => ({ __isRef: true, blockId: "b1", name });

/** Every input the drawer offers, filled in, so only a setting can be the reason a throw happens. */
const configured: BlockData = {
  ...initBlockData(),
  abundanceRef: ref("abundance"),
  filterRef: ref("filter"),
  concentrationColumnRef: ref("concentration"),
  binColumnRef: ref("bin"),
  binFluorescenceColumnRef: ref("binFluorescence"),
  antigenColumnRef: ref("antigen"),
  sortFractionColumnRef: ref("sortFraction"),
  targetAntigen: "CD19",
  customBlockLabel: "CD19 titration",
};

const hillConfigured: BlockData = { ...configured, equation: "hill_gate_index" };

/**
 * The exact bytes `deriveArgs` projected before the plateau-floor gate and the `FitSettingsDraft`
 * split. The args are the block's cache key: a change in these strings re-runs every project that
 * already holds results, so this pair has to stay fixed unless a change to the fit intends it.
 */
const LOG_FLUORESCENCE_ARGS_JSON =
  '{"abundanceRef":{"__isRef":true,"blockId":"b1","name":"abundance"},"concentrationColumnRef":{"__isRef":true,"blockId":"b1","name":"concentration"},"binColumnRef":{"__isRef":true,"blockId":"b1","name":"bin"},"antigenColumnRef":{"__isRef":true,"blockId":"b1","name":"antigen"},"sortFractionColumnRef":{"__isRef":true,"blockId":"b1","name":"sortFraction"},"targetAntigen":"CD19","binFluorescenceColumnRef":{"__isRef":true,"blockId":"b1","name":"binFluorescence"},"equation":"log_fluorescence","plateauFloor":0.667,"minReadsPerConcentration":3,"minConcentrationPoints":5,"r2GateMin":0.8,"kdSeGateMax":1,"traceLabel":"CD19 titration","inputFilter":"{\\"__isRef\\":true,\\"blockId\\":\\"b1\\",\\"name\\":\\"filter\\"}"}';

const HILL_GATE_INDEX_ARGS_JSON =
  '{"abundanceRef":{"__isRef":true,"blockId":"b1","name":"abundance"},"concentrationColumnRef":{"__isRef":true,"blockId":"b1","name":"concentration"},"binColumnRef":{"__isRef":true,"blockId":"b1","name":"bin"},"antigenColumnRef":{"__isRef":true,"blockId":"b1","name":"antigen"},"sortFractionColumnRef":{"__isRef":true,"blockId":"b1","name":"sortFraction"},"targetAntigen":"CD19","minReadsPerConcentration":3,"minConcentrationPoints":5,"r2ThresholdGood":0.8,"r2ThresholdFailed":0.5,"nMin":0.5,"nMax":2,"hookEffectThresholdBin":0.2,"hookEffectThresholdNoBin":0.02,"hookEffectMinReads":20,"traceLabel":"CD19 titration","inputFilter":"{\\"__isRef\\":true,\\"blockId\\":\\"b1\\",\\"name\\":\\"filter\\"}"}';

describe("a required numeric setting the scientist cleared disables Run", () => {
  // One row per setting the selected equation reads. A throw out of `deriveArgs` is what the
  // platform reads as "args invalid", which is what greys the Run button out.
  const cases: { setting: keyof BlockData; data: BlockData; message: string }[] = [
    { setting: "plateauFloor", data: configured, message: "Plateau floor is required" },
    {
      setting: "minReadsPerConcentration",
      data: hillConfigured,
      message: "Min reads per concentration is required",
    },
    {
      setting: "minConcentrationPoints",
      data: hillConfigured,
      message: "Min concentration points is required",
    },
    {
      setting: "r2ThresholdGood",
      data: hillConfigured,
      message: "R² threshold (Good) is required",
    },
    {
      setting: "r2ThresholdFailed",
      data: hillConfigured,
      message: "R² threshold (Failed) is required",
    },
    { setting: "nMin", data: hillConfigured, message: "Hill coefficient nMin is required" },
    { setting: "nMax", data: hillConfigured, message: "Hill coefficient nMax is required" },
    {
      setting: "hookEffectThresholdBin",
      data: hillConfigured,
      message: "Hook effect signal-drop threshold (bin mode) is required",
    },
    {
      setting: "hookEffectMinReads",
      data: hillConfigured,
      message: "Min reads for hook check is required",
    },
  ];

  it.each(cases)("$setting", ({ setting, data, message }) => {
    expect(() => deriveArgs({ ...data, [setting]: undefined })).toThrow(message);
  });

  it("hookEffectThresholdNoBin, in the frequency mode that reads it", () => {
    expect(() =>
      deriveArgs({
        ...hillConfigured,
        binColumnRef: undefined,
        sortFractionColumnRef: undefined,
        hookEffectThresholdNoBin: undefined,
      }),
    ).toThrow("Hook effect signal-drop threshold (frequency mode) is required");
  });
});

describe("a setting the selected equation never reads cannot disable Run", () => {
  it("log fluorescence runs with every hill_gate_index setting cleared", () => {
    const cleared: BlockData = {
      ...configured,
      r2ThresholdGood: undefined,
      r2ThresholdFailed: undefined,
      nMin: undefined,
      nMax: undefined,
      hookEffectThresholdBin: undefined,
      hookEffectThresholdNoBin: undefined,
      hookEffectMinReads: undefined,
    };
    expect(JSON.stringify(deriveArgs(cleared))).toBe(LOG_FLUORESCENCE_ARGS_JSON);
  });

  it("hill_gate_index runs with the plateau floor cleared", () => {
    expect(JSON.stringify(deriveArgs({ ...hillConfigured, plateauFloor: undefined }))).toBe(
      HILL_GATE_INDEX_ARGS_JSON,
    );
  });

  it("hill_gate_index runs with the inactive mode's hook threshold cleared", () => {
    // Bin mode, so the frequency-mode threshold is hidden in the drawer. The args still carry it,
    // at its default, because the workflow's schema has the field.
    expect(() =>
      deriveArgs({ ...hillConfigured, hookEffectThresholdNoBin: undefined }),
    ).not.toThrow();
  });
});

describe("valid settings project the bytes they always projected", () => {
  it("under the log fluorescence fit", () => {
    expect(JSON.stringify(deriveArgs(configured))).toBe(LOG_FLUORESCENCE_ARGS_JSON);
  });

  it("under the mean bin fit", () => {
    expect(JSON.stringify(deriveArgs(hillConfigured))).toBe(HILL_GATE_INDEX_ARGS_JSON);
  });
});

describe("the plateau floor's own range", () => {
  // 0 is the one value `min-value="0"` lets through that the fit cannot use: a floor of 0 asks
  // the curve to rise by nothing.
  it.each([0, -0.1, 1.5])("rejects %s", (plateauFloor) => {
    expect(() => deriveArgs({ ...configured, plateauFloor })).toThrow(
      "Plateau floor must be above 0 and at most 1",
    );
  });

  it.each([Number.MIN_VALUE, 0.4, 1])("accepts %s", (plateauFloor) => {
    expect(() => deriveArgs({ ...configured, plateauFloor })).not.toThrow();
  });
});

describe("the trace label", () => {
  it("is the scientist's own label when set", () => {
    expect(deriveArgs(configured).traceLabel).toBe("CD19 titration");
  });

  it("names the block and the target antigen otherwise", () => {
    expect(deriveArgs({ ...configured, customBlockLabel: "" }).traceLabel).toBe(
      "Tite-Seq Analysis (CD19)",
    );
  });

  it("names the block alone without a target antigen", () => {
    const data = { ...configured, antigenColumnRef: undefined, targetAntigen: undefined };
    expect(deriveArgs({ ...data, customBlockLabel: "" }).traceLabel).toBe("Tite-Seq Analysis");
  });

  it("does not depend on the subtitle's pool-derived part", () => {
    const data = { ...configured, customBlockLabel: "" };
    expect(deriveArgs({ ...data, subtitleLabel: "Abundance" })).toEqual(deriveArgs(data));
  });
});
