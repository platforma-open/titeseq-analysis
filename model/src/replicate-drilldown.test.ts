import type { PColumnSpec, PlRef } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { BlockData } from "./index";
import { dataModel, deriveArgs, initBlockData, replicateDrilldownFilter } from "./index";

const ref = (name: string): PlRef => ({ __isRef: true, blockId: "b1", name });

const VARIANT_AXIS = {
  name: "pl7.app/vdj/clonotypeKey",
  type: "String",
  domain: { "pl7.app/vdj/clonotypingRunId": "run1" },
  annotations: { "pl7.app/label": "Clonotype key" },
} as const;
const REPLICATE_AXIS = {
  name: "pl7.app/vdj/replicate",
  type: "String",
  domain: { "pl7.app/block": "blk" },
} as const;

const perReplicate = (name: string): PColumnSpec => ({
  kind: "PColumn",
  name,
  valueType: "Double",
  domain: { "pl7.app/block": "blk" },
  axesSpec: [VARIANT_AXIS, REPLICATE_AXIS],
});

const frame = [
  perReplicate("pl7.app/vdj/kd"),
  perReplicate("pl7.app/vdj/curveFitR2"),
  perReplicate("pl7.app/vdj/kdReplicateState"),
];

describe("the replicate drill-down filter", () => {
  it("narrows the frame to the clicked variant on its own axis, domain included", () => {
    expect(replicateDrilldownFilter(frame, ["v-17"])).toEqual({
      type: "and",
      filters: [
        {
          type: "patternEquals",
          column: {
            type: "axis",
            id: {
              name: "pl7.app/vdj/clonotypeKey",
              type: "String",
              domain: { "pl7.app/vdj/clonotypingRunId": "run1" },
            },
          },
          value: "v-17",
        },
      ],
    });
  });

  it("finds the variant axis whichever position it holds", () => {
    const swapped = frame.map((s) => ({ ...s, axesSpec: [REPLICATE_AXIS, VARIANT_AXIS] }));
    const filter = replicateDrilldownFilter(swapped, ["v-17"]);
    expect(JSON.stringify(filter)).toContain("pl7.app/vdj/clonotypeKey");
    expect(JSON.stringify(filter)).not.toContain("pl7.app/vdj/replicate");
  });

  it("builds nothing while no variant is selected", () => {
    expect(replicateDrilldownFilter(frame, undefined)).toBeUndefined();
    expect(replicateDrilldownFilter(frame, [])).toBeUndefined();
  });

  it("builds nothing without a per-replicate Kd,app column", () => {
    expect(replicateDrilldownFilter(frame.slice(1), ["v-17"])).toBeUndefined();
  });
});

describe("opening the drill-down never stales the block", () => {
  const configured: BlockData = {
    ...initBlockData(),
    abundanceRef: ref("abundance"),
    concentrationColumnRef: ref("concentration"),
    binColumnRef: ref("bin"),
    replicateColumnRef: ref("replicate"),
  };

  it("the selected key and its grid state stay out of args", () => {
    const opened: BlockData = {
      ...configured,
      replicateDrilldownKey: ["v-17"],
      replicateDrilldownTableState: initBlockData().tableState,
    };
    expect(JSON.stringify(deriveArgs(opened))).toBe(JSON.stringify(deriveArgs(configured)));
  });

  it("a fresh block starts with the drill-down closed", () => {
    expect(initBlockData().replicateDrilldownKey).toBeUndefined();
  });

  it("a stored v4 project without the fields reads as closed", () => {
    expect("replicateDrilldownKey" in configured).toBe(false);
    const migrated = dataModel.migrate({ version: "v4", data: configured });
    expect(migrated.version).toBe("v7");
    expect(migrated.data.replicateDrilldownKey).toBeUndefined();
  });
});
