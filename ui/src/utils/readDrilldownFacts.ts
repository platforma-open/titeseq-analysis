import { drilldownColumnKey } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { getRawPlatformaInstance, isLabelColumn, pTableValue } from "@platforma-sdk/model";
import type { PTableHandle } from "@platforma-sdk/model";

export type DrilldownFacts = {
  /** The variant's label; undefined when the table has no label column or no row. */
  label: string | undefined;
  /** Keys (`drilldownColumnKey`) of the columns with no value in any row. */
  emptyColumns: string[];
};

/**
 * What the drill-down panel reads from its table: the first row's label, and the columns that
 * hold no value. A drill-down table holds one variant, so row 0 names it. The model cannot read
 * either: the clonotype axis value is an opaque key, and the columns are Parquet-partitioned.
 */
export async function readDrilldownFacts(handle: PTableHandle): Promise<DrilldownFacts> {
  const driver = getRawPlatformaInstance().pFrameDriver;
  const specs = await driver.getSpec(handle);
  const { rows } = await driver.getShape(handle);
  if (rows === 0) return { label: undefined, emptyColumns: [] };
  const columns = specs.flatMap((s, i) => (s.type === "column" ? [{ spec: s.spec, i }] : []));
  const data = await driver.getData(
    handle,
    columns.map((c) => c.i),
    { offset: 0, length: rows },
  );
  let label: string | undefined;
  const emptyColumns: string[] = [];
  columns.forEach(({ spec }, at) => {
    const values = Array.from({ length: rows }, (_, row) => pTableValue(data[at], row));
    if (isLabelColumn(spec)) {
      if (typeof values[0] === "string") label = values[0];
    } else if (values.every((v) => v === null || v === "")) {
      emptyColumns.push(drilldownColumnKey(spec));
    }
  });
  return { label, emptyColumns };
}
