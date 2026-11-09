import { getRawPlatformaInstance, isLabelColumn, pTableValue } from "@platforma-sdk/model";
import type { PTableHandle } from "@platforma-sdk/model";

/**
 * The first row's value in a table's label column. A drill-down table holds one variant, so row 0
 * names it. The model cannot resolve the label: the clonotype axis value is an opaque key, and its
 * label column is Parquet-partitioned.
 */
export async function readTableLabel(handle: PTableHandle): Promise<string | undefined> {
  const driver = getRawPlatformaInstance().pFrameDriver;
  const specs = await driver.getSpec(handle);
  const iLabel = specs.findIndex((s) => s.type === "column" && isLabelColumn(s.spec));
  if (iLabel < 0) return undefined;
  const { rows } = await driver.getShape(handle);
  if (rows === 0) return undefined;
  const [labels] = await driver.getData(handle, [iLabel], { offset: 0, length: 1 });
  const label = pTableValue(labels, 0);
  return typeof label === "string" ? label : undefined;
}
