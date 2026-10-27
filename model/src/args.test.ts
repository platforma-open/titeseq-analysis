import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockData } from "./index";
import { deriveArgs, initBlockData, routeFitModel } from "./index";

const ref = (name: string): PlRef => ({ __isRef: true, blockId: "b1", name });

const binMode = (overrides: Partial<BlockData> = {}): BlockData => ({
  ...initBlockData(),
  abundanceRef: ref("abundance"),
  concentrationColumnRef: ref("concentration"),
  binColumnRef: ref("bin"),
  cellCountColumnRef: ref("cells"),
  ...overrides,
});

// Settings the UI hides under the Adams model, each set to a value its own check rejects.
const invalidHillSettings: [string, Partial<BlockData>][] = [
  ["nMin >= nMax", { nMin: 3, nMax: 1 }],
  ["Failed R² above Good", { r2ThresholdGood: 0.4, r2ThresholdFailed: 0.6 }],
  ["R² out of range", { r2ThresholdGood: 1.5 }],
];

describe("routing the curve model on the bound columns", () => {
  it("measures when the cells sorted into each sample are known", () => {
    expect(routeFitModel(binMode())).toBe("poisson_log");
  });

  it("measures when sort fractions stand in for cell counts", () => {
    expect(
      routeFitModel(
        binMode({ cellCountColumnRef: undefined, sortFractionColumnRef: ref("fractions") }),
      ),
    ).toBe("poisson_log");
  });

  it("ranks when bins are bound but the cells are unknown", () => {
    expect(routeFitModel(binMode({ cellCountColumnRef: undefined }))).toBe("hill_free_n");
  });

  it("ranks in frequency mode, where the linear Hill is the only model", () => {
    expect(routeFitModel(binMode({ binColumnRef: undefined, cellCountColumnRef: undefined }))).toBe(
      "hill_free_n",
    );
  });

  it("lets the override win over the routing", () => {
    expect(deriveArgs(binMode({ fitModelOverride: "hill_free_n" })).fitModel).toBe("hill_free_n");
    expect(
      deriveArgs(binMode({ cellCountColumnRef: undefined, fitModelOverride: "poisson_log" }))
        .fitModel,
    ).toBe("poisson_log");
  });
});

describe("the args projection under the Adams model", () => {
  it.each(invalidHillSettings)("ignores a hidden Hill-only setting: %s", (_, bad) => {
    // A user binds cell counts, so the field they left invalid under the Hill fit is now hidden.
    expect(() => deriveArgs(binMode(bad))).not.toThrow();
  });

  it("sends the same args whatever the hidden Hill-only settings hold", () => {
    // Editing a hidden setting must not mark the block stale.
    const plain = deriveArgs(binMode());
    const edited = deriveArgs(binMode({ nMin: 0.9, r2ThresholdGood: 0.95 }));
    expect(edited).toEqual(plain);
  });
});

describe("the args projection when the linear Hill model runs", () => {
  it.each(invalidHillSettings)("rejects an invalid Hill setting in bin mode: %s", (_, bad) => {
    expect(() => deriveArgs(binMode({ cellCountColumnRef: undefined, ...bad }))).toThrow();
  });

  it("rejects nMin >= nMax in frequency mode, where Hill is the only model", () => {
    const data = binMode({
      binColumnRef: undefined,
      cellCountColumnRef: undefined,
      nMin: 3,
      nMax: 1,
    });
    expect(() => deriveArgs(data)).toThrow("nMin must be < nMax");
  });

  it("carries the user's Hill settings through", () => {
    const args = deriveArgs(binMode({ cellCountColumnRef: undefined, nMin: 0.7, nMax: 1.8 }));
    expect([args.nMin, args.nMax]).toEqual([0.7, 1.8]);
  });
});

describe("the args projection and the inputs that no longer block a run", () => {
  it("runs with bins and no cells at all", () => {
    expect(() =>
      deriveArgs(binMode({ cellCountColumnRef: undefined, sortFractionColumnRef: undefined })),
    ).not.toThrow();
  });

  it("carries the optional FACS refinements through", () => {
    const args = deriveArgs(
      binMode({
        binFluorescenceColumnRef: ref("meanFluor"),
        gateWidthColumnRef: ref("gateWidth"),
      }),
    );
    expect(args.binFluorescenceColumnRef).toEqual(ref("meanFluor"));
    expect(args.gateWidthColumnRef).toEqual(ref("gateWidth"));
  });

  it("refuses a FACS refinement without a bin column", () => {
    expect(() =>
      deriveArgs(
        binMode({
          binColumnRef: undefined,
          cellCountColumnRef: undefined,
          binFluorescenceColumnRef: ref("meanFluor"),
        }),
      ),
    ).toThrow("Measured bin fluorescence requires a FACS bin column");
  });

  it("rejects an amplitude prior outside (0, 1)", () => {
    expect(() => deriveArgs(binMode({ amplitudePriorQ: 1 }))).toThrow("Amplitude prior");
    expect(() => deriveArgs(binMode({ amplitudePriorQ: 0.65 }))).not.toThrow();
  });

  it("defaults the hook detector off and the Hill exponent pinned", () => {
    const args = deriveArgs(binMode());
    expect(args.hookDetectorEnabled).toBe(false);
    expect(args.pinHillExponent).toBe(true);
  });
});
