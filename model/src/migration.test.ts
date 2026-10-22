import { describe, expect, it } from "vitest";
import type { BlockDataV1 } from "./index";
import { initBlockData, migrateV1ToV2 } from "./index";

describe("the v1 -> v2 migration", () => {
  // A v1 project: everything a block had before the Adams 2016 model, configured.
  const {
    cellCountColumnRef,
    fitModel,
    binWidthDecades,
    minPlateauFraction,
    baselineFluorescence,
    allowUncorrected,
    kdSdGood,
    ...v1
  } = { ...initBlockData(), minReadsPerConcentration: 7, customBlockLabel: "run 4" };
  void [
    cellCountColumnRef,
    fitModel,
    binWidthDecades,
    minPlateauFraction,
    baselineFluorescence,
    allowUncorrected,
    kdSdGood,
  ];

  it("keeps an existing project on the linear Hill it ran with", () => {
    // It has no cell counts, which the new default model requires; a rerun must not start failing.
    expect(migrateV1ToV2(v1 as BlockDataV1).fitModel).toBe("hill_free_n");
  });

  it("keeps every v1 setting and fills the new ones with the block defaults", () => {
    const v2 = migrateV1ToV2(v1 as BlockDataV1);
    expect(v2.minReadsPerConcentration).toBe(7);
    expect(v2.customBlockLabel).toBe("run 4");
    const fresh = initBlockData();
    expect(v2.binWidthDecades).toBe(fresh.binWidthDecades);
    expect(v2.minPlateauFraction).toBe(fresh.minPlateauFraction);
    expect(v2.kdSdGood).toBe(fresh.kdSdGood);
    expect(v2.allowUncorrected).toBe(false);
    expect(v2.cellCountColumnRef).toBeUndefined();
  });

  it("starts a new block on the Adams 2016 model", () => {
    expect(initBlockData().fitModel).toBe("poisson_log");
  });
});
