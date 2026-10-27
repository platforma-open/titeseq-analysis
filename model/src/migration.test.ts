import { describe, expect, it } from "vitest";
import type { BlockDataV1, BlockDataV2 } from "./index";
import { initBlockData, migrateV1ToV2, migrateV2ToV3, routeFitModel } from "./index";

// A v2 project: what the block stored once the Adams 2016 model existed and was chosen by hand.
const v2Project = (fitModel: "poisson_log" | "hill_free_n"): BlockDataV2 => {
  const {
    fitModelOverride,
    binFluorescenceColumnRef,
    gateWidthColumnRef,
    amplitudePriorQ,
    hookDetectorEnabled,
    pinHillExponent,
    ...rest
  } = { ...initBlockData(), minReadsPerConcentration: 7, customBlockLabel: "run 4" };
  void [
    fitModelOverride,
    binFluorescenceColumnRef,
    gateWidthColumnRef,
    amplitudePriorQ,
    hookDetectorEnabled,
    pinHillExponent,
  ];
  return { ...rest, fitModel, allowUncorrected: false };
};

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
  } = v2Project("poisson_log");
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
    expect(v2.cellCountColumnRef).toBeUndefined();
  });
});

describe("the v2 -> v3 migration", () => {
  // `migrate` fast-paths a project already at the latest version, so a project that passed
  // through v1 -> v2 would keep its pinned model forever if v3 did not drop it.
  it("drops the stored model so the project is routed like a new one", () => {
    const v3 = migrateV2ToV3(v2Project("hill_free_n"));
    expect(v3.fitModelOverride).toBeUndefined();
    expect("fitModel" in v3).toBe(false);
  });

  it("leaves a migrated project routed on what it actually has bound", () => {
    const v3 = migrateV2ToV3(v2Project("hill_free_n"));
    // No bin column was configured in the fixture, so routing gives the ranking fit.
    expect(routeFitModel(v3)).toBe("hill_free_n");
  });

  it("drops the uncorrected-fit flag, which no longer gates anything", () => {
    expect("allowUncorrected" in migrateV2ToV3(v2Project("poisson_log"))).toBe(false);
  });

  it("fills the new settings with the block defaults", () => {
    const v3 = migrateV2ToV3(v2Project("poisson_log"));
    expect(v3.hookDetectorEnabled).toBe(false);
    expect(v3.pinHillExponent).toBe(true);
    expect(v3.amplitudePriorQ).toBeUndefined();
  });

  it("keeps every setting the project already had", () => {
    const v3 = migrateV2ToV3(v2Project("poisson_log"));
    expect(v3.minReadsPerConcentration).toBe(7);
    expect(v3.customBlockLabel).toBe("run 4");
  });
});
