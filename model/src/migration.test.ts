import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockArgs, BlockData, BlockDataV1 } from "./index";
import {
  HILL_GATE_INDEX_DEFAULTS,
  dataModel,
  deriveArgs,
  initBlockData,
  migrateV1ToV2,
  migrateV2ToV3,
} from "./index";

const ref = (name: string): PlRef => ({ __isRef: true, blockId: "b1", name });

/**
 * A project stored by the block as published — v2.4.0's `BlockData`, which had no equation to
 * name. The settings are deliberately off their defaults: a migration that quietly reseeds them
 * would pass against a project that happened to sit at every default.
 */
const v240Project: BlockDataV1 = {
  abundanceRef: ref("abundance"),
  filterRef: ref("filter"),
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
  defaultBlockLabel: "CD19 - Abundance",
  customBlockLabel: "CD19 titration",
  tableState: initBlockData().tableState,
  graphStateTitrationCurves: initBlockData().graphStateTitrationCurves,
  graphStateKDHistogram: initBlockData().graphStateKDHistogram,
  graphStateAffinityVsFit: initBlockData().graphStateAffinityVsFit,
  settingsOpen: true,
};

/** The project as the block will hold it after the upgrade, through the real version chain. */
const upgraded = (): BlockData => {
  const migrated = dataModel.migrate({ version: "v1", data: v240Project });
  expect(migrated.version).toBe("v3");
  return migrated.data;
};

/** The args the model projects from a given state — what the staleness gate compares. */
const argsOf = (data: BlockData): BlockArgs => deriveArgs(data);

describe("a project created on v2.4.0 or earlier", () => {
  it("is pinned to the equation it ran, by name", () => {
    const data = upgraded();
    // Not `toBe` alone: an absent field would read `undefined`, and `undefined` defaults forward
    // to the new equation the moment anything reads it.
    expect(Object.hasOwn(data, "equation")).toBe(true);
    expect(data.equation).toBe("hill_gate_index");
  });

  it("does not pick up the equation a new block gets", () => {
    expect(initBlockData().equation).toBe("log_fluorescence");
    expect(upgraded().equation).not.toBe(initBlockData().equation);
  });

  it("keeps every setting, reference and piece of view state it stored", () => {
    const data = upgraded();
    for (const [field, value] of Object.entries(v240Project)) {
      expect(data[field as keyof BlockDataV1]).toEqual(value);
    }
  });

  it("projects the args v2.4.0 projected, so the upgrade does not stale it", () => {
    // v2.4.0's args lambda is not in the tree, but its output is reconstructible: the v2 data
    // model adds `equation` and `plateauFloor` to args and nothing else, so the bytes must match.
    const args = argsOf(upgraded());
    expect(Object.hasOwn(args, "equation")).toBe(false);
    expect(Object.hasOwn(args, "plateauFloor")).toBe(false);
  });

  it("carries a plateau floor its equation never reads", () => {
    expect(upgraded().plateauFloor).toBe(2 / 3);
  });

  it("carries an R² gate its equation never reads, and its args name none", () => {
    expect(upgraded()).toMatchObject({ r2GateEnabled: false, r2GateMin: 0.8 });
    expect(Object.hasOwn(argsOf(upgraded()), "r2GateMin")).toBe(false);
  });
});

describe("a project created on a v2 development build", () => {
  const v2Project = migrateV1ToV2(v240Project);

  it("arrives at v3 with the R² gate off, so it rates its fits as it did", () => {
    const migrated = dataModel.migrate({
      version: "v2",
      data: { ...v2Project, equation: "log_fluorescence" },
    });
    expect(migrated.version).toBe("v3");
    expect(migrated.data).toMatchObject({ r2GateEnabled: false, r2GateMin: 0.8 });
  });

  it("moves nothing else", () => {
    const { r2GateEnabled: _e, r2GateMin: _m, ...rest } = migrateV2ToV3(v2Project);
    expect(rest).toEqual(v2Project);
  });
});

describe("a block created on the v2 data model", () => {
  it("fits the log-fluorescence equation", () => {
    expect(initBlockData().equation).toBe("log_fluorescence");
    expect(initBlockData().plateauFloor).toBe(2 / 3);
  });

  it("sends its equation and floor by name, so nothing has to be inferred", () => {
    const args = argsOf({
      ...initBlockData(),
      abundanceRef: ref("abundance"),
      concentrationColumnRef: ref("concentration"),
      // The default equation fits gate fluorescence, so a new block needs gates.
      binColumnRef: ref("bin"),
    });
    expect(args).toMatchObject({ equation: "log_fluorescence", plateauFloor: 2 / 3 });
  });

  it("cannot run the default equation on frequency-only data", () => {
    // No-bin mode survives only on v2.4.0's equation. Run stays disabled with a reason
    // rather than the job failing in Python after submission.
    expect(() =>
      argsOf({
        ...initBlockData(),
        abundanceRef: ref("abundance"),
        concentrationColumnRef: ref("concentration"),
        binColumnRef: undefined,
      }),
    ).toThrow("requires a FACS bin column");
  });

  it("runs v2.4.0's equation on frequency-only data, as it always did", () => {
    expect(() =>
      argsOf({
        ...initBlockData(),
        equation: "hill_gate_index",
        abundanceRef: ref("abundance"),
        concentrationColumnRef: ref("concentration"),
        binColumnRef: undefined,
      }),
    ).not.toThrow();
  });
});

describe("the args carry only the settings the selected equation reads", () => {
  const configured: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
    binFluorescenceColumnRef: ref("binFluorescence"),
  };

  it("log fluorescence omits every hill_gate_index setting", () => {
    const args = argsOf(configured);
    for (const field of Object.keys(HILL_GATE_INDEX_DEFAULTS)) {
      expect(Object.hasOwn(args, field), field).toBe(false);
    }
    expect(args).toMatchObject({ binFluorescenceColumnRef: ref("binFluorescence") });
  });

  it("log fluorescence runs with an invalid hill_gate_index setting the drawer hides", () => {
    expect(() => argsOf({ ...configured, nMin: 3, nMax: 2 })).not.toThrow();
  });

  it("a hidden hill_gate_index setting does not change log fluorescence args", () => {
    expect(argsOf({ ...configured, r2ThresholdGood: 0.95 })).toEqual(argsOf(configured));
  });

  it("hill_gate_index omits the plateau floor, the R² gate and gate fluorescence", () => {
    const args = argsOf({ ...configured, equation: "hill_gate_index", r2GateEnabled: true });
    for (const field of [
      "plateauFloor",
      "r2GateEnabled",
      "r2GateMin",
      "binFluorescenceColumnRef",
    ]) {
      expect(Object.hasOwn(args, field), field).toBe(false);
    }
  });
});

describe("the R² gate", () => {
  const configured: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
  };

  it("is off for a new block, and its args name no gate", () => {
    expect(configured.r2GateEnabled).toBe(false);
    expect(Object.hasOwn(argsOf(configured), "r2GateMin")).toBe(false);
  });

  it("an off gate's threshold does not change the args", () => {
    expect(argsOf({ ...configured, r2GateMin: 0.5 })).toEqual(argsOf(configured));
  });

  it("sends its threshold while on", () => {
    expect(argsOf({ ...configured, r2GateEnabled: true })).toMatchObject({ r2GateMin: 0.8 });
  });

  it("disables Run while on with its threshold cleared or out of range", () => {
    expect(() => argsOf({ ...configured, r2GateEnabled: true, r2GateMin: undefined })).toThrow(
      "R² gate threshold is required",
    );
    expect(() => argsOf({ ...configured, r2GateEnabled: true, r2GateMin: 1.2 })).toThrow(
      "between 0 and 1",
    );
  });

  it("an off gate with its threshold cleared still runs", () => {
    expect(() => argsOf({ ...configured, r2GateMin: undefined })).not.toThrow();
  });
});

describe("the migration function on its own", () => {
  it("writes the equation into the data, not only into the args", () => {
    expect(migrateV1ToV2(v240Project).equation).toBe("hill_gate_index");
  });

  it("is a no-op for everything else", () => {
    const { equation: _e, plateauFloor: _p, ...rest } = migrateV1ToV2(v240Project);
    expect(rest).toEqual(v240Project);
  });
});
