import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockData } from "./index";
import { binRanks, defaultBinOrder, deriveArgs, initBlockData, unrankedBins } from "./index";

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

  it("keeps a skipped gate number's gap under the default order", () => {
    // Gate 3 has no samples: gate 4 stays gate 4, as a block without an order reads it.
    expect(binRanks(defaultBinOrder(["4", "1", "2"]))).toEqual({ "1": 1, "2": 2, "4": 4 });
    expect(binRanks(["2", "3", "4"])).toEqual({ "2": 2, "3": 3, "4": 4 });
  });

  it("ranks by position once the numbers are reordered or are not positive integers", () => {
    expect(binRanks(["4", "2", "1"])).toEqual({ "4": 1, "2": 2, "1": 3 });
    expect(binRanks(["0", "1", "2"])).toEqual({ "0": 1, "1": 2, "2": 3 });
    expect(binRanks(["1", "2.5"])).toEqual({ "1": 1, "2.5": 2 });
  });

  it("disables Run on a label column with no order", () => {
    expect(() => deriveArgs({ ...configured, binColumnIsText: true })).toThrow(
      "Set the FACS bin order",
    );
    expect(
      deriveArgs({ ...configured, binColumnIsText: true, binOrder: ["LOW", "HIGH"] }),
    ).toHaveProperty("binRanks", { LOW: 1, HIGH: 2 });
    // An integer column runs without an order.
    expect(deriveArgs({ ...configured, binColumnIsText: false })).not.toHaveProperty("binRanks");
  });

  it("names the column's values the order does not rank", () => {
    expect(unrankedBins(["1", "2"], ["1", "2", "3"])).toEqual(["3"]);
    expect(unrankedBins(["1", "2", "3"], ["1", "2"])).toEqual([]);
    expect(unrankedBins(undefined, ["1"])).toEqual([]);
    expect(unrankedBins(["1"], undefined)).toEqual([]);
  });

  it("seeds numbers numerically and labels alphabetically", () => {
    expect(defaultBinOrder(["10", "2", "1", "2"])).toEqual(["1", "2", "10"]);
    expect(defaultBinOrder(["mid", "high", "low"])).toEqual(["high", "low", "mid"]);
    expect(binRanks(["b", "a"])).toEqual({ b: 1, a: 2 });
  });
});

describe("the synonymous-mode input queries", () => {
  it("are left out only when the variant axis is known not to be nucleotide", () => {
    expect(deriveArgs({ ...configured, variantAxisIsNucleotide: false })).toHaveProperty(
      "synonymousInputs",
      false,
    );
    expect(deriveArgs({ ...configured, variantAxisIsNucleotide: true })).not.toHaveProperty(
      "synonymousInputs",
    );
    // A block that stored no snapshot keeps its args.
    expect(deriveArgs(configured)).not.toHaveProperty("synonymousInputs");
  });
});
