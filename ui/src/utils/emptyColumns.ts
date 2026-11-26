import { REPLICATE_AXIS } from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  canonicalizeJson,
  createPlDataTableStateV2,
  getPTableColumnId,
  getRawPlatformaInstance,
  isColumnHidden,
  isColumnOptional,
  pTableValue,
} from "@platforma-sdk/model";
import type {
  PTableColumnId,
  PTableColumnSpec,
  PTableHandle,
  PlDataTableColumnsMeta,
  PlDataTableStateV2Normalized,
} from "@platforma-sdk/model";

/** A column's visibility as the grid resolves it: the model's rule first, then the spec. */
function visibilityOf(
  spec: PTableColumnSpec,
  meta: PlDataTableColumnsMeta | undefined,
): "default" | "optional" | "hidden" {
  if (spec.type === "axis") return "default";
  const rule = meta?.columns[spec.id]?.visibility;
  if (rule !== undefined) return rule;
  if (isColumnHidden(spec.spec)) return "hidden";
  if (isColumnOptional(spec.spec)) return "optional";
  return "default";
}

const isEmptyValue = (value: unknown): boolean => value === null || value === "";

/**
 * The columns a variant's replicate table starts with hidden: the model's optional columns, and
 * every column it would show that holds no value in any replicate. `values[i]` holds column i's
 * cells; axes are never hidden.
 */
export function defaultHiddenColumns(
  specs: PTableColumnSpec[],
  values: unknown[][],
  meta: PlDataTableColumnsMeta | undefined,
): PTableColumnId[] {
  return specs.flatMap((spec, i) => {
    if (spec.type === "axis") return [];
    const visibility = visibilityOf(spec, meta);
    if (visibility === "hidden") return [];
    if (visibility === "optional" || values[i].every(isEmptyValue)) {
      return [getPTableColumnId(spec)];
    }
    return [];
  });
}

/**
 * A fresh table state whose entry for `sourceId` starts with `hidden` columns hidden. Everything
 * else is the state the grid builds for a source it has not seen.
 */
export function seededTableState(
  sourceId: string,
  hidden: PTableColumnId[],
): PlDataTableStateV2Normalized {
  const base = createPlDataTableStateV2();
  const hiddenColIds = hidden.map((id) => canonicalizeJson<PTableColumnId>(id));
  return {
    ...base,
    stateCache: [
      {
        sourceId,
        gridState: { columnVisibility: { hiddenColIds } },
        sorting: null,
        hiddenColIds: hidden,
        sheetsState: [],
        filtersState: null,
        defaultFiltersState: null,
      },
    ],
    // The model builds the table from the cache entry; the params name the source only.
    pTableParams: { ...base.pTableParams, sourceId },
  };
}

/**
 * Reads a variant's replicate table and returns the columns it starts with hidden. Undefined when
 * the table holds another variant: the output still carries the previous selection's table until
 * the model recomputes it.
 */
export async function readDefaultHiddenColumns(
  handle: PTableHandle,
  variantKey: string,
  meta: PlDataTableColumnsMeta | undefined,
): Promise<PTableColumnId[] | undefined> {
  const driver = getRawPlatformaInstance().pFrameDriver;
  const specs = await driver.getSpec(handle);
  const { rows } = await driver.getShape(handle);
  if (rows === 0) return [];
  const vectors = await driver.getData(
    handle,
    specs.map((_, i) => i),
    { offset: 0, length: rows },
  );
  const values = vectors.map((vector) =>
    Array.from({ length: rows }, (_, row) => pTableValue(vector, row, { na: null })),
  );
  const iVariant = specs.findIndex((s) => s.type === "axis" && s.spec.name !== REPLICATE_AXIS);
  if (iVariant < 0) {
    throw new Error("The replicate table has no variant axis besides the replicate axis.");
  }
  if (String(values[iVariant][0]) !== variantKey) return undefined;
  return defaultHiddenColumns(specs, values, meta);
}
