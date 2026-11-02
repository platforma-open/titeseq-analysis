import type { PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockArgs, BlockData, BlockDataV1, BlockDataV2, BlockDataV3 } from "./index";
import {
  HILL_GATE_INDEX_DEFAULTS,
  dataModel,
  deriveArgs,
  initBlockData,
  migrateV1ToV2,
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
  expect(migrated.version).toBe("v4");
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

  it("hill_gate_index omits the plateau floor and gate fluorescence", () => {
    const args = argsOf({ ...configured, equation: "hill_gate_index" });
    expect(Object.hasOwn(args, "plateauFloor")).toBe(false);
    expect(Object.hasOwn(args, "binFluorescenceColumnRef")).toBe(false);
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

describe("a project stored on the v2 data model", () => {
  const {
    ntTableState: _n,
    replicateTableState: _r,
    graphStateReplicates: _g,
    ...v2Project
  }: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
    customBlockLabel: "CD19 titration",
    settingsOpen: true,
  };
  const stored: BlockDataV2 = v2Project;
  const migrated = () => dataModel.migrate({ version: "v2", data: stored });

  it("gains an empty nucleotide variant table state", () => {
    expect(migrated().version).toBe("v4");
    expect(migrated().data.ntTableState).toEqual(initBlockData().ntTableState);
  });

  it("keeps everything it stored", () => {
    const {
      ntTableState: _added,
      replicateTableState: _addedTable,
      graphStateReplicates: _addedGraph,
      ...rest
    } = migrated().data;
    expect(rest).toEqual(stored);
  });

  it("projects the args it projected", () => {
    expect(JSON.stringify(argsOf(migrated().data))).toBe(
      JSON.stringify(argsOf({ ...initBlockData(), ...stored })),
    );
  });
});

describe("a project stored on the v3 data model", () => {
  const {
    replicateTableState: _r,
    graphStateReplicates: _g,
    ...v3Project
  }: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
    customBlockLabel: "CD19 titration",
    settingsOpen: true,
  };
  const stored: BlockDataV3 = v3Project;
  const migrated = () => dataModel.migrate({ version: "v3", data: stored });

  it("gains the per-replicate view state and nothing else", () => {
    expect(migrated().version).toBe("v4");
    expect(migrated().data.replicateTableState).toEqual(initBlockData().replicateTableState);
    expect(migrated().data.graphStateReplicates).toEqual(initBlockData().graphStateReplicates);
  });

  it("keeps everything it stored", () => {
    const { replicateTableState: _t, graphStateReplicates: _s, ...rest } = migrated().data;
    expect(rest).toEqual(stored);
  });

  it("projects not one byte of replicate argument", () => {
    const args = argsOf(migrated().data);
    expect(Object.keys(args).filter((k) => k.startsWith("replicate"))).toEqual([]);
    expect(JSON.stringify(args)).toBe(JSON.stringify(argsOf({ ...initBlockData(), ...stored })));
  });
});

describe("the pages", () => {
  const hrefs = (synonymous: boolean, replicates = false) =>
    sectionsFor(synonymous, replicates).map((s) => s.href);

  it("show the nucleotide variant table in synonymous mode, after the main table", () => {
    expect(hrefs(true)).toEqual([
      "/",
      "/nt-variants",
      "/titration-curves",
      "/kd-distribution",
      "/affinity-vs-fit",
    ]);
  });

  it("leave every other input with the four pages it had", () => {
    expect(sectionsFor(false)).toEqual([
      { type: "link", href: "/", label: "Table" },
      { type: "link", href: "/titration-curves", label: "Titration Curves" },
      { type: "link", href: "/kd-distribution", label: "Kd Distribution" },
      { type: "link", href: "/affinity-vs-fit", label: "Affinity vs Fit Quality" },
    ]);
  });

  it("show the per-replicate page only when replicates were combined", () => {
    expect(hrefs(false, true)).toEqual([
      "/",
      "/replicates",
      "/titration-curves",
      "/kd-distribution",
      "/affinity-vs-fit",
    ]);
    expect(hrefs(false, false)).not.toContain("/replicates");
  });

  it("put the per-replicate page after the nucleotide table in synonymous mode", () => {
    expect(hrefs(true, true)).toEqual([
      "/",
      "/nt-variants",
      "/replicates",
      "/titration-curves",
      "/kd-distribution",
      "/affinity-vs-fit",
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
      replicateErrorModel: "phillips" as const,
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
      replicateErrorModel: "sd_over_sqrt_n",
      replicateCensoredMajority: 0.5,
    });
  });

  it("leaves both filters out of the args until they are set", () => {
    const args = argsOf({ ...configured, replicateColumnRef: ref("replicate") });
    expect(Object.hasOwn(args, "replicateFitErrorCutoff")).toBe(false);
    expect(Object.hasOwn(args, "replicateR2Floor")).toBe(false);
  });

  it("sends each filter once it is set", () => {
    const args = argsOf({
      ...configured,
      replicateColumnRef: ref("replicate"),
      replicateFitErrorCutoff: 1,
      replicateR2Floor: 0.8,
    });
    expect(args).toMatchObject({ replicateFitErrorCutoff: 1, replicateR2Floor: 0.8 });
  });

  it("refuses a filter or a majority the software would reject", () => {
    const withColumn = { ...configured, replicateColumnRef: ref("replicate") };
    expect(() => argsOf({ ...withColumn, replicateFitErrorCutoff: 0 })).toThrow();
    expect(() => argsOf({ ...withColumn, replicateR2Floor: 1.5 })).toThrow();
    expect(() => argsOf({ ...withColumn, replicateCensoredMajority: 0.4 })).toThrow();
  });

  it("is unset on a new block and on an upgraded project", () => {
    expect(initBlockData().replicateColumnRef).toBeUndefined();
    expect(upgraded().replicateColumnRef).toBeUndefined();
    expect(replicateKeys(upgraded())).toEqual([]);
  });
});
