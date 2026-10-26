import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockData } from "./index";
import { deriveArgs, initBlockData } from "./index";

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

describe("the args projection under the Adams model", () => {
  it.each(invalidHillSettings)("ignores a hidden Hill-only setting: %s", (_, bad) => {
    // A user switches from linear Hill to Adams; the field they left invalid is now hidden.
    expect(() => deriveArgs(binMode({ fitModel: "poisson_log", ...bad }))).not.toThrow();
  });

  it("sends the same args whatever the hidden Hill-only settings hold", () => {
    // Editing a hidden setting must not mark the block stale.
    const plain = deriveArgs(binMode({ fitModel: "poisson_log" }));
    const edited = deriveArgs(
      binMode({ fitModel: "poisson_log", nMin: 0.9, r2ThresholdGood: 0.95 }),
    );
    expect(edited).toEqual(plain);
  });
});

describe("the args projection when the linear Hill model runs", () => {
  it.each(invalidHillSettings)("rejects an invalid Hill setting in bin mode: %s", (_, bad) => {
    expect(() => deriveArgs(binMode({ fitModel: "hill_free_n", ...bad }))).toThrow();
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
    const args = deriveArgs(binMode({ fitModel: "hill_free_n", nMin: 0.7, nMax: 1.8 }));
    expect([args.nMin, args.nMax]).toEqual([0.7, 1.8]);
  });
});
