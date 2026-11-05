import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockData } from "./index";
import { binRanks, defaultBinOrder, deriveArgs, initBlockData } from "./index";

const ref = (name: string): PlRef => ({ __isRef: true, blockId: "b1", name });

const configured: BlockData = {
  ...initBlockData(),
  abundanceRef: ref("abundance"),
  concentrationColumnRef: ref("concentration"),
  binColumnRef: ref("bin"),
};

describe("the FACS bin order", () => {
  it("projects no ranks when no order was set, so older blocks keep their args", () => {
    expect(deriveArgs(configured)).not.toHaveProperty("binRanks");
  });

  it("ranks each kept value by its position, from 1", () => {
    const args = deriveArgs({ ...configured, binOrder: ["LOW", "MID", "HIGH"] });
    expect(args).toHaveProperty("binRanks", { LOW: 1, MID: 2, HIGH: 3 });
  });

  it("projects no ranks without a bin column", () => {
    // The mean bin fit runs without bins; the log fluorescence fit needs them.
    const args = deriveArgs({
      ...configured,
      equation: "hill_gate_index",
      binColumnRef: undefined,
      binOrder: ["1", "2"],
    });
    expect(args).not.toHaveProperty("binRanks");
  });

  it("disables Run when every bin was removed", () => {
    expect(() => deriveArgs({ ...configured, binOrder: [] })).toThrow("Keep at least one FACS bin");
  });

  it("seeds numbers numerically and labels alphabetically", () => {
    expect(defaultBinOrder(["10", "2", "1", "2"])).toEqual(["1", "2", "10"]);
    expect(defaultBinOrder(["mid", "high", "low"])).toEqual(["high", "low", "mid"]);
    expect(binRanks(["b", "a"])).toEqual({ b: 1, a: 2 });
  });
});
