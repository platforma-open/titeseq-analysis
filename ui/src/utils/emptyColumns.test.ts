import type {
  PColumnSpec,
  PObjectId,
  PTableColumnSpec,
  PlDataTableColumnsMeta,
} from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import { defaultHiddenColumns, seededTableState } from "./emptyColumns";

const variantAxis = { name: "pl7.app/variant", type: "String" as const };

function axis(): PTableColumnSpec {
  return { type: "axis", id: variantAxis, spec: variantAxis };
}

function column(
  id: string,
  annotations: Record<string, string> = {},
): Extract<PTableColumnSpec, { type: "column" }> {
  const spec: PColumnSpec = {
    kind: "PColumn",
    name: `pl7.app/${id}`,
    valueType: "Double",
    axesSpec: [variantAxis],
    annotations,
  };
  return { type: "column", id: id as PObjectId, spec };
}

const columnId = (id: string) => ({ type: "column", id });

describe("defaultHiddenColumns", () => {
  it("hides a column with no value in any replicate", () => {
    const specs = [axis(), column("kd"), column("reason")];
    const values = [
      ["v1", "v1"],
      [1.5, null],
      [null, null],
    ];
    expect(defaultHiddenColumns(specs, values, undefined)).toEqual([columnId("reason")]);
  });

  it("counts an empty string as no value", () => {
    const specs = [axis(), column("reason")];
    expect(defaultHiddenColumns(specs, [["v1"], [""]], undefined)).toEqual([columnId("reason")]);
  });

  it("keeps the model's optional columns hidden even with values", () => {
    const specs = [axis(), column("se", { "pl7.app/table/visibility": "optional" })];
    expect(defaultHiddenColumns(specs, [["v1"], [0.2]], undefined)).toEqual([columnId("se")]);
  });

  it("follows the model's rule over the spec annotation", () => {
    const se = column("se", { "pl7.app/table/visibility": "optional" });
    const specs = [axis(), se];
    const meta: PlDataTableColumnsMeta = {
      columns: { [se.id]: { visibility: "default" as const } },
      axes: {},
    };
    expect(defaultHiddenColumns(specs, [["v1"], [0.2]], meta)).toEqual([]);
  });

  it("leaves out columns the grid never shows, empty or not", () => {
    const specs = [axis(), column("label", { "pl7.app/table/visibility": "hidden" })];
    expect(defaultHiddenColumns(specs, [["v1"], [null]], undefined)).toEqual([]);
  });

  it("never hides an axis", () => {
    expect(defaultHiddenColumns([axis()], [[null]], undefined)).toEqual([]);
  });
});

describe("seededTableState", () => {
  it("holds one entry for the source with its hidden columns; the params name the source", () => {
    const hidden = [{ type: "column" as const, id: "reason" as PObjectId }];
    const state = seededTableState('"v1"', hidden);
    expect(state.version).toBe(9);
    expect(state.stateCache).toEqual([
      {
        sourceId: '"v1"',
        gridState: { columnVisibility: { hiddenColIds: ['{"id":"reason","type":"column"}'] } },
        sorting: null,
        hiddenColIds: hidden,
        sheetsState: [],
        filtersState: null,
        defaultFiltersState: null,
      },
    ]);
    expect(state.pTableParams).toEqual({
      sourceId: '"v1"',
      hiddenColIds: null,
      sorting: null,
      filters: null,
      defaultFilters: null,
    });
  });
});
