/**
 * The variants behind the row the variant panel opened. Read from the graph frame in the UI, so
 * opening a row writes nothing to the block's data.
 */
import {
  CONCENTRATION_AXIS,
  REPLICATE_AXIS,
  VARIANT_CHART,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import type {
  AxisSpec,
  PColumnIdAndSpec,
  PFrameHandle,
  PObjectId,
  PTableColumnId,
} from "@platforma-sdk/model";
import { getAxisId, getRawPlatformaInstance } from "@platforma-sdk/model";

/** Each fitted nucleotide variant's amino-acid variant: the panel finds a row's synonymous variants by it. */
export const AA_VARIANT_KEY = "pl7.app/vdj/aaVariantKey";
const LABEL = "pl7.app/label";
const FITTED_SIGNAL = VARIANT_CHART.fitted;
const BLOCK_DOMAIN = "pl7.app/block";

/** One variant on the curves' axis: its key and its label. */
export type PanelVariant = { key: string; label: string };

const sameAxis = (a: AxisSpec, b: AxisSpec) =>
  JSON.stringify(getAxisId(a)) === JSON.stringify(getAxisId(b));

/**
 * The column `name` that sits on `axis` alone. A graph frame also holds the result pool, with the
 * columns of other blocks of this kind: `block`, when given, keeps this block's own.
 */
function onAxis(cols: PColumnIdAndSpec[], name: string, axis: AxisSpec, block?: string) {
  return cols.find(
    (c) =>
      c.spec.name === name &&
      c.spec.axesSpec.length === 1 &&
      sameAxis(c.spec.axesSpec[0], axis) &&
      (block === undefined || c.spec.domain?.[BLOCK_DOMAIN] === block),
  );
}

/** A frame output's handle while the output is ready; undefined while computing or failed. */
export const handleOf = <T>(output: { ok: true; value: T } | { ok: false } | undefined) =>
  output?.ok ? output.value : undefined;

/** The rows of `columns` where `filter` equals `key`: the key axis values, then each column's. */
async function rowsWhere(
  frame: PFrameHandle,
  columns: PObjectId[],
  filter: PTableColumnId,
  key: string,
): Promise<{ keys: string[]; values: Map<string, unknown[]> }> {
  const [primary, ...rest] = columns;
  // The specs come from reactive model outputs, and a reactive proxy cannot cross to the driver.
  const plainFilter: PTableColumnId = JSON.parse(JSON.stringify(filter));
  const response = await getRawPlatformaInstance().pFrameDriver.calculateTableData(frame, {
    src: {
      type: "outer",
      primary: { type: "column", column: primary },
      secondary: rest.map((column) => ({ type: "column" as const, column })),
    },
    filters: [
      {
        type: "bySingleColumnV2",
        column: plainFilter,
        predicate: { operator: "Equal", reference: key },
      },
    ],
    sorting: [],
  });
  const keys: string[] = [];
  const values = new Map<string, unknown[]>();
  for (const column of response) {
    const data = Array.from(column.data.data as ArrayLike<unknown>);
    if (column.spec.type === "axis") {
      if (keys.length === 0) keys.push(...data.map(String));
    } else {
      values.set(String(column.spec.id), data);
    }
  }
  return { keys, values };
}

/** The curves' variant axis and the block that wrote them; undefined without the fitted signal. */
async function curveAxis(frame: PFrameHandle) {
  const cols = await getRawPlatformaInstance().pFrameDriver.listColumns(frame);
  const signal = cols.find((c) => c.spec.name === FITTED_SIGNAL);
  const axis = signal?.spec.axesSpec.find(
    (a) => a.name !== REPLICATE_AXIS && a.name !== CONCENTRATION_AXIS,
  );
  const block = signal?.spec.domain?.[BLOCK_DOMAIN];
  return axis === undefined || block === undefined ? undefined : { cols, axis, block };
}

/**
 * Synonymous mode: the synonymous variants of the opened amino-acid variant that drew a curve, on
 * the curves' own axis, by label. Empty on a run that did not link them.
 */
export async function loadSynonymousVariants(
  frame: PFrameHandle,
  rowKey: string,
): Promise<PanelVariant[]> {
  const found = await curveAxis(frame);
  if (found === undefined) return [];
  const { cols, axis, block } = found;
  const aaKey = onAxis(cols, AA_VARIANT_KEY, axis, block);
  if (aaKey === undefined) return [];
  const label = onAxis(cols, LABEL, axis);
  const { keys, values } = await rowsWhere(
    frame,
    [aaKey.columnId, ...(label ? [label.columnId] : [])],
    { type: "column", id: aaKey.columnId },
    rowKey,
  );
  return keys
    .map((key, i) => {
      const text = label && values.get(label.columnId)?.[i];
      return { key, label: typeof text === "string" && text !== "" ? text : key };
    })
    .sort((a, b) => a.label.localeCompare(b.label));
}

/**
 * The label of one variant on the curves' axis: the result pool's label column, which the graph
 * frame carries. Undefined where the input has no label column or the variant no label.
 */
export async function loadCurveVariantLabel(
  frame: PFrameHandle,
  key: string,
): Promise<string | undefined> {
  const found = await curveAxis(frame);
  if (found === undefined) return undefined;
  const label = onAxis(found.cols, LABEL, found.axis);
  if (label === undefined) return undefined;
  const { values } = await rowsWhere(
    frame,
    [label.columnId],
    { type: "axis", id: getAxisId(found.axis) },
    key,
  );
  const text = values.get(label.columnId)?.[0];
  return typeof text === "string" && text !== "" ? text : undefined;
}

/**
 * Whether the variant `key` has any point of the curve column `signalName`, among the frame's
 * columns `cols`. False without the column.
 */
export async function hasCurve(
  frame: PFrameHandle,
  cols: PColumnIdAndSpec[],
  signalName: string,
  key: string,
): Promise<boolean> {
  const signal = cols.find((c) => c.spec.name === signalName);
  const axis = signal?.spec.axesSpec.find(
    (a) => a.name !== REPLICATE_AXIS && a.name !== CONCENTRATION_AXIS,
  );
  if (signal === undefined || axis === undefined) return false;
  const { keys } = await rowsWhere(
    frame,
    [signal.columnId],
    { type: "axis", id: getAxisId(axis) },
    key,
  );
  return keys.length > 0;
}
