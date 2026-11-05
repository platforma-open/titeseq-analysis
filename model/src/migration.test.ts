import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockArgs, BlockData, BlockDataV1, BlockDataV3, BlockDataV5 } from "./index";
import {
  HILL_GATE_INDEX_DEFAULTS,
  dataModel,
  deriveArgs,
  initBlockData,
  migrateV1ToV2,
  migrateV2ToV3,
  migrateV3ToV4,
  migrateV4ToV5,
  migrateV5ToV6,
  sectionsFor,
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
  expect(migrated.version).toBe("v7");
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
    expect(upgraded().plateauFloor).toBe(0.667);
  });

  it("arrives with no R² limit, and its args name none", () => {
    expect(upgraded().r2GateMin).toBeUndefined();
    expect(Object.hasOwn(upgraded(), "r2GateEnabled")).toBe(false);
    expect(Object.hasOwn(argsOf(upgraded()), "r2GateMin")).toBe(false);
  });
});

describe("a project created on a v2 development build", () => {
  const v2Project = migrateV1ToV2(v240Project);

  it("arrives with no R² limit, so it rates its fits as it did", () => {
    const migrated = dataModel.migrate({
      version: "v2",
      data: { ...v2Project, equation: "log_fluorescence" },
    });
    expect(migrated.version).toBe("v7");
    expect(migrated.data.r2GateMin).toBeUndefined();
  });

  it("moves nothing else", () => {
    const { r2GateEnabled: _e, r2GateMin: _m, ...rest } = migrateV2ToV3(v2Project);
    expect(rest).toEqual(v2Project);
  });
});

describe("a block created on the v2 data model", () => {
  it("fits the log-fluorescence equation", () => {
    expect(initBlockData().equation).toBe("log_fluorescence");
    expect(initBlockData().plateauFloor).toBe(0.667);
  });

  it("sends its equation and floor by name, so nothing has to be inferred", () => {
    const args = argsOf({
      ...initBlockData(),
      abundanceRef: ref("abundance"),
      concentrationColumnRef: ref("concentration"),
      // The default equation fits gate fluorescence, so a new block needs gates.
      binColumnRef: ref("bin"),
    });
    expect(args).toMatchObject({ equation: "log_fluorescence", plateauFloor: 0.667 });
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
    const args = argsOf({ ...configured, equation: "hill_gate_index", r2GateMin: 0.8 });
    for (const field of ["plateauFloor", "r2GateMin", "binFluorescenceColumnRef"]) {
      expect(Object.hasOwn(args, field), field).toBe(false);
    }
  });
});

describe("the R² and SE limits", () => {
  const configured: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
  };

  it("start a new block at Minimum R² 0.8 and Maximum Kd,app SE 1", () => {
    expect(argsOf(configured)).toMatchObject({ r2GateMin: 0.8, kdSeGateMax: 1 });
  });

  it("send no R² limit when it is cleared", () => {
    expect(Object.hasOwn(argsOf({ ...configured, r2GateMin: undefined }), "r2GateMin")).toBe(false);
  });

  it("send each limit that is set, on its own", () => {
    expect(argsOf({ ...configured, r2GateMin: 0.6 })).toMatchObject({ r2GateMin: 0.6 });
    const seOnly = argsOf({ ...configured, r2GateMin: undefined, kdSeGateMax: 1 });
    expect(seOnly).toMatchObject({ kdSeGateMax: 1 });
    expect(Object.hasOwn(seOnly, "r2GateMin")).toBe(false);
  });

  it("disable Run on a limit out of range", () => {
    expect(() => argsOf({ ...configured, r2GateMin: 1.2 })).toThrow("between 0 and 1");
    expect(() => argsOf({ ...configured, r2GateMin: -0.1 })).toThrow("between 0 and 1");
    expect(() => argsOf({ ...configured, kdSeGateMax: 0 })).toThrow(
      "Maximum Kd,app SE must be above 0",
    );
  });

  it("never sends the complete-titration rule, which is withheld", () => {
    for (const requireAllConcentrations of [undefined, false, true]) {
      const args = argsOf({ ...configured, requireAllConcentrations });
      expect(Object.hasOwn(args, "requireAllConcentrations")).toBe(false);
    }
  });
});

describe("the synonymous weighting", () => {
  const configured: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
  };
  const equations = ["log_fluorescence", "hill_gate_index"] as const;

  it.each(equations)("left unset or none, %s projects the args it always did", (equation) => {
    const base = { ...configured, equation };
    expect(Object.hasOwn(argsOf(base), "synonymousWeighting")).toBe(false);
    expect(JSON.stringify(argsOf({ ...base, synonymousWeighting: "none" }))).toBe(
      JSON.stringify(argsOf(base)),
    );
  });

  it.each(equations)("set to reads, %s sends it last", (equation) => {
    const args = argsOf({ ...configured, equation, synonymousWeighting: "reads" });
    expect(args.synonymousWeighting).toBe("reads");
    expect(Object.keys(args).at(-1)).toBe("synonymousWeighting");
  });

  it("is unset on a new block and on an upgraded project", () => {
    expect(initBlockData().synonymousWeighting).toBeUndefined();
    expect(upgraded().synonymousWeighting).toBeUndefined();
  });
});

describe("a project stored on the v6 data model", () => {
  const { ntTableState: _n, ...v6Project }: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
    customBlockLabel: "CD19 titration",
    settingsOpen: true,
  };
  // A v6 project 2.5.0 wrote with a current error-model name.
  const stored = v6Project;
  const migrated = () => dataModel.migrate({ version: "v6", data: stored });

  it("renames the error model 2.5.0 stored as phillips", () => {
    const legacy = dataModel.migrate({
      version: "v6",
      data: { ...stored, replicateErrorModel: "phillips" },
    });
    expect(legacy.data.replicateErrorModel).toBe("pop_sd_over_sqrt_n_plus_1");
    const other = dataModel.migrate({
      version: "v6",
      data: { ...stored, replicateErrorModel: "sd_over_sqrt_n" },
    });
    expect(other.data.replicateErrorModel).toBe("sd_over_sqrt_n");
  });

  it("gains an empty nucleotide variant table state", () => {
    expect(migrated().version).toBe("v7");
    expect(migrated().data.ntTableState).toEqual(initBlockData().ntTableState);
  });

  it("keeps everything it stored", () => {
    const { ntTableState: _added, ...rest } = migrated().data;
    expect(rest).toEqual(stored);
  });

  it("projects the args it projected", () => {
    expect(JSON.stringify(argsOf(migrated().data))).toBe(
      JSON.stringify(argsOf({ ...stored, ntTableState: initBlockData().ntTableState })),
    );
  });
});

describe("the pages", () => {
  const hrefs = (synonymous: boolean) => sectionsFor({ synonymous }).map((s) => s.href);

  it("show the nucleotide variant table in synonymous mode, after the main table", () => {
    expect(hrefs(true)).toEqual([
      "/",
      "/nt-variants",
      "/synonymous-noise",
      "/titration-edges",
      "/read-depth",
      "/titration-curves",
      "/kd-distribution",
      "/affinity-vs-fit",
    ]);
  });

  it("show the parent, synonymous and replicate pages only where they apply", () => {
    expect(
      sectionsFor({ synonymous: true, parentMode: true, replicates: true }).map((s) => s.href),
    ).toEqual([
      "/",
      "/nt-variants",
      "/vs-parent",
      "/synonymous-noise",
      "/replicate-noise",
      "/titration-edges",
      "/read-depth",
      "/titration-curves",
      "/kd-distribution",
      "/affinity-vs-fit",
    ]);
  });

  it("give every other input the summary pages beside the four it had", () => {
    expect(sectionsFor({ synonymous: false })).toEqual([
      { type: "link", href: "/", label: "Main" },
      { type: "link", href: "/titration-edges", label: "Titration Edges" },
      { type: "link", href: "/read-depth", label: "Read Depth" },
      { type: "link", href: "/titration-curves", label: "Titration Curves" },
      { type: "link", href: "/kd-distribution", label: "Kd Distribution" },
      { type: "link", href: "/affinity-vs-fit", label: "Affinity vs Fit Quality" },
    ]);
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

describe("a project stored on the v3 data model", () => {
  const stored: BlockDataV3 = migrateV2ToV3(migrateV1ToV2(v240Project));
  const migrated = () => dataModel.migrate({ version: "v3", data: stored });

  it("gains unset replicates and the per-replicate view state", () => {
    expect(migrated().version).toBe("v7");
    expect(migrated().data).toMatchObject({
      replicateErrorModel: "pop_sd_over_sqrt_n_plus_1",
      replicateCensoredMajority: 0.5,
    });
    expect(migrated().data.replicateColumnRef).toBeUndefined();
    expect(migrated().data.replicateTableState).toEqual(initBlockData().replicateTableState);
    expect(migrated().data.graphStateReplicates).toEqual(initBlockData().graphStateReplicates);
  });

  it("moves nothing else", () => {
    const {
      replicateErrorModel: _m,
      replicateCensoredMajority: _c,
      replicateTableState: _t,
      graphStateReplicates: _g,
      ...rest
    } = migrateV3ToV4(stored);
    expect(rest).toEqual(stored);
  });

  it("projects not one byte of replicate argument", () => {
    const args = argsOf(migrated().data);
    expect(Object.keys(args).filter((k) => k.startsWith("replicate"))).toEqual([]);
  });
});

describe("a project stored on the v4 data model", () => {
  const stored = (minConcentrationPoints: number | undefined): BlockDataV5 => ({
    ...migrateV3ToV4(migrateV2ToV3(migrateV1ToV2(v240Project))),
    minConcentrationPoints,
  });
  const migrated = (n: number | undefined) => dataModel.migrate({ version: "v4", data: stored(n) });

  it.each([3, 4])("raises a point count of %i to 5", (n) => {
    expect(migrated(n).version).toBe("v7");
    expect(migrated(n).data.minConcentrationPoints).toBe(5);
  });

  it.each([5, 8])("keeps a point count of %i", (n) => {
    expect(migrated(n).data.minConcentrationPoints).toBe(n);
  });

  it("keeps a cleared point count cleared", () => {
    expect(migrated(undefined).data.minConcentrationPoints).toBeUndefined();
  });

  it("moves nothing else", () => {
    const { minConcentrationPoints: _n, ...rest } = migrateV4ToV5(stored(3));
    const { minConcentrationPoints: _o, ...before } = stored(3);
    expect(rest).toEqual(before);
  });
});

describe("the sections", () => {
  it("show only Main before the first run", () => {
    expect(sectionsFor(undefined)).toEqual([{ type: "link", href: "/", label: "Main" }]);
  });

  it("add only the noise page when replicates were combined", () => {
    // Per-replicate values are read by clicking a variant in the main table; the noise page is
    // the one replicate page.
    expect(sectionsFor({ synonymous: false, replicates: true }).map((s) => s.href)).toEqual([
      "/",
      "/replicate-noise",
      "/titration-edges",
      "/read-depth",
      "/titration-curves",
      "/kd-distribution",
      "/affinity-vs-fit",
    ]);
  });
});

describe("the sample replicates", () => {
  const configured: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
  };
  const equations = ["log_fluorescence", "hill_gate_index"] as const;
  const replicateKeys = (data: BlockData) =>
    Object.keys(argsOf(data)).filter((k) => k.startsWith("replicate"));

  it.each(equations)("%s sends no replicate byte without a column", (equation) => {
    expect(replicateKeys({ ...configured, equation })).toEqual([]);
  });

  it("a stored setting alone changes nothing: the column is what turns the mode on", () => {
    const tuned = {
      ...configured,
      replicateErrorModel: "sd_over_sqrt_n" as const,
      replicateFitErrorCutoff: 1,
      replicateR2Floor: 0.8,
      replicateCensoredMajority: 0.9,
    };
    expect(JSON.stringify(argsOf(tuned))).toBe(JSON.stringify(argsOf(configured)));
  });

  it("the column brings the error model and the censoring majority with it", () => {
    const args = argsOf({ ...configured, replicateColumnRef: ref("replicate") });
    expect(args).toMatchObject({
      replicateColumnRef: ref("replicate"),
      replicateErrorModel: "pop_sd_over_sqrt_n_plus_1",
      replicateCensoredMajority: 0.5,
    });
  });

  it("leaves both filters out of the args until they are set", () => {
    const args = argsOf({ ...configured, replicateColumnRef: ref("replicate") });
    expect(Object.hasOwn(args, "replicateFitErrorCutoff")).toBe(false);
    expect(Object.hasOwn(args, "replicateR2Floor")).toBe(false);
  });

  it("pins every hidden replicate setting to its default while the section is hidden", () => {
    // A block that stored these before the section went away must not run them unseen.
    const stored = {
      ...configured,
      replicateColumnRef: ref("replicate"),
      replicateErrorModel: "sd_over_sqrt_n" as const,
      replicateFitErrorCutoff: 1,
      replicateR2Floor: 0.8,
      replicateCensoredMajority: 0.9,
    };
    expect(JSON.stringify(argsOf(stored))).toBe(
      JSON.stringify(argsOf({ ...configured, replicateColumnRef: ref("replicate") })),
    );
  });

  it("does not block Run on a hidden setting the operator cannot see", () => {
    const withColumn = { ...configured, replicateColumnRef: ref("replicate") };
    expect(() => argsOf({ ...withColumn, replicateFitErrorCutoff: 0 })).not.toThrow();
    expect(() => argsOf({ ...withColumn, replicateR2Floor: 1.5 })).not.toThrow();
    expect(() => argsOf({ ...withColumn, replicateCensoredMajority: undefined })).not.toThrow();
  });

  it("is unset on a new block and on an upgraded project", () => {
    expect(initBlockData().replicateColumnRef).toBeUndefined();
    expect(upgraded().replicateColumnRef).toBeUndefined();
    expect(replicateKeys(upgraded())).toEqual([]);
  });
});

describe("a project stored on the v5 data model", () => {
  const base = (): BlockDataV5 => ({
    ...migrateV4ToV5(migrateV3ToV4(migrateV2ToV3(migrateV1ToV2(v240Project)))),
    equation: "log_fluorescence",
    r2GateMin: 0.85,
    kdSeGateMax: 1.5,
  });
  const migrated = (r2GateEnabled: boolean) =>
    dataModel.migrate({ version: "v5", data: { ...base(), r2GateEnabled } });
  it("with the checkbox off, clears both limits", () => {
    expect(migrated(false).version).toBe("v7");
    expect(migrated(false).data.r2GateMin).toBeUndefined();
    expect(migrated(false).data.kdSeGateMax).toBeUndefined();
  });

  it("with the checkbox on, keeps both limits", () => {
    expect(migrated(true).data).toMatchObject({ r2GateMin: 0.85, kdSeGateMax: 1.5 });
  });

  it("projects the args it projected, so it does not re-run", () => {
    // v5 sent both limits with the checkbox on, and neither with it off.
    const off = argsOf(migrated(false).data);
    expect(Object.hasOwn(off, "r2GateMin")).toBe(false);
    expect(Object.hasOwn(off, "kdSeGateMax")).toBe(false);
    expect(argsOf(migrated(true).data)).toEqual({ ...off, r2GateMin: 0.85, kdSeGateMax: 1.5 });
  });

  it("drops the checkbox and moves nothing else", () => {
    const { r2GateEnabled: _e, ...rest } = base();
    const data = migrateV5ToV6({ ...base(), r2GateEnabled: true });
    expect(Object.hasOwn(data, "r2GateEnabled")).toBe(false);
    expect(data).toEqual(rest);
  });
});
