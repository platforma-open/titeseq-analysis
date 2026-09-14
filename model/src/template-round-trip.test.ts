import { kind } from "@platforma-open/platforma-open.titeseq-analysis.kind";
import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockData } from "./index";
import { deriveTemplateParams, initBlockData } from "./index";

/**
 * Export a block's data as a template would, then create a block from it.
 *
 * The `JSON` hop is deliberate: a template is a file, so anything that survives only in memory
 * is not actually carried. What comes back is a fresh block's data, which is what a scientist
 * applying the template gets.
 */
const roundTrip = (data: BlockData): BlockData =>
  initBlockData(
    kind.parseInitializationParams(JSON.parse(JSON.stringify(deriveTemplateParams(data)))),
  );

const ref = (name: string): PlRef => ({ __isRef: true, blockId: "b1", name });

/** A fully configured block: every field the contract carries, none of them at its default. */
const configured: BlockData = {
  ...initBlockData(),
  abundanceRef: ref("abundance"),
  concentrationColumnRef: ref("concentration"),
  binColumnRef: ref("bin"),
  antigenColumnRef: ref("antigen"),
  sortFractionColumnRef: ref("sortFraction"),
  targetAntigen: "CD19",
  minReadsPerConcentration: 7,
  minConcentrationPoints: 8,
  r2ThresholdGood: 0.9,
  r2ThresholdFailed: 0.3,
  nMin: 0.75,
  nMax: 3,
  hookEffectThresholdBin: 0.35,
  hookEffectThresholdNoBin: 0.05,
  hookEffectMinReads: 50,
  customBlockLabel: "CD19 titration",
};

describe("the template round trip", () => {
  it("carries every field the contract names", () => {
    const restored = roundTrip(configured);
    for (const field of Object.keys(deriveTemplateParams(configured)) as (keyof BlockData)[]) {
      expect(restored[field]).toEqual(configured[field]);
    }
  });

  it("is idempotent — a second pass changes nothing", () => {
    expect(roundTrip(roundTrip(configured))).toEqual(roundTrip(configured));
  });

  it("carries the values a `??` default would swallow", () => {
    const restored = roundTrip({
      ...configured,
      r2ThresholdFailed: 0,
      hookEffectThresholdNoBin: 0,
      customBlockLabel: "",
    });
    expect(restored.r2ThresholdFailed).toBe(0);
    expect(restored.hookEffectThresholdNoBin).toBe(0);
    expect(restored.customBlockLabel).toBe("");
  });

  it("carries a half-configured block — an abundance column and nothing else", () => {
    const restored = roundTrip({ ...initBlockData(), abundanceRef: ref("abundance") });
    expect(restored.abundanceRef).toEqual(ref("abundance"));
    expect(restored.concentrationColumnRef).toBeUndefined();
    expect(restored.minConcentrationPoints).toBe(5);
  });

  it("gives an untouched block back unchanged", () => {
    expect(roundTrip(initBlockData())).toEqual(initBlockData());
  });

  it("does not carry the derived label, the drawer flag or view state", () => {
    const restored = roundTrip({
      ...configured,
      defaultBlockLabel: "CD19 - Abundance",
      settingsOpen: true,
    });
    // Re-derived from the antigen and abundance labels by a `watchEffect` on mount.
    expect(restored.defaultBlockLabel).toBe("Tite-Seq Analysis");
    expect(restored.settingsOpen).toBe(false);
    expect(restored.tableState).toEqual(initBlockData().tableState);
  });
});
