import { assertParamsObject, defineBlockKind } from "@platforma-sdk/block-kind";
import type { PlRef } from "@platforma-sdk/model";
import { isPlRef } from "@platforma-sdk/model";
import { isString } from "es-toolkit";
import { isNumber } from "es-toolkit/compat";
import { name, version } from "../package.json" with { type: "json" };

/**
 * This block's init-params contract — what a creator or a project template supplies to seed a
 * new instance: the columns a titration is read from, the fit's tuning, and the user's own label.
 *
 * Excluded on purpose:
 *   * `defaultBlockLabel` — derived from the target antigen and the abundance column's label, and
 *     rewritten by a `watchEffect` the moment the UI mounts. A template value for it would be
 *     overwritten before anyone could read it. Only the scientist's own `customBlockLabel` is
 *     worth restoring.
 *   * `settingsOpen` — whether the settings drawer is showing, which `onMounted` opens by itself
 *     for a block with no abundance column picked.
 *   * `tableState` and the three graph states — view state.
 *
 * Every field is optional: a block may be created without a template, and a template need not
 * set all of them.
 */
export type BlockParams = {
  abundanceRef?: PlRef;
  concentrationColumnRef?: PlRef;
  binColumnRef?: PlRef;
  antigenColumnRef?: PlRef;
  sortFractionColumnRef?: PlRef;
  targetAntigen?: string;
  minReadsPerConcentration?: number;
  minConcentrationPoints?: number;
  r2ThresholdGood?: number;
  r2ThresholdFailed?: number;
  nMin?: number;
  nMax?: number;
  hookEffectThresholdBin?: number;
  hookEffectThresholdNoBin?: number;
  hookEffectMinReads?: number;
  customBlockLabel?: string;
};

/**
 * The contract at runtime, for params arriving from a template file rather than typed code. An
 * absent field is always allowed — every param is optional and the block's own default takes
 * over — so each guard runs only on what is present. Keys the contract does not name are dropped
 * by never being read.
 *
 * Rules that span two fields are not checked here: that `nMin` is below `nMax`, that the failed
 * R² threshold is at or below the good one, that a target antigen is named once an antigen column
 * is picked, that a FACS sort fraction needs a bin column. Every one of those is a state the
 * settings drawer leaves behind mid-edit, and `.args()` already refuses them when the block tries
 * to run. Refusing them here would mean a block the scientist can reach by hand cannot be carried
 * by a template.
 */
function parseInitializationParams(value: unknown): BlockParams {
  assertParamsObject(value);

  const params: Record<string, unknown> = {};
  for (const [field, { is, must }] of Object.entries(CONTRACT)) {
    const v = value[field];
    if (v === undefined) continue;
    if (!is(v)) throw new Error(`'${field}' must be ${must}.`);
    params[field] = v;
  }
  return params as BlockParams;
}

// Identity (`name`/`version`) comes from this package's own `package.json`, so the on-wire
// `{name}@{version}` reference can never drift from what npm publishes; the bundler inlines the
// JSON import.
export const kind = defineBlockKind<BlockParams>({
  name,
  version,
  parseInitializationParams,
});

// ---------------------------------------------------------------------------
// Internals
// ---------------------------------------------------------------------------

type Guard<T> = (v: unknown) => v is T;

/** A guard plus how to finish the sentence "'field' must be …". */
type Check<T> = { is: Guard<T>; must: string };

function check<T>(is: Guard<T>, must: string): Check<T> {
  return { is, must };
}

/** A real number. `isNumber` alone admits `NaN` and the infinities, which no setting means. */
const isFiniteNumber: Guard<number> = (v): v is number => isNumber(v) && Number.isFinite(v);

/** A number at or above `min`, the way the args lambda words the same bound. */
function isNumberFrom(min: number): Guard<number> {
  return (v): v is number => isFiniteNumber(v) && v >= min;
}

const isFraction: Guard<number> = (v): v is number => isFiniteNumber(v) && v >= 0 && v <= 1;

/**
 * The runtime half of the contract. The `satisfies` clause is what stops it drifting: every
 * field `BlockParams` declares must appear here, and each guard must narrow to that field's own
 * type — so adding a param without a check stops compiling.
 *
 * The bounds mirror the args lambda exactly, so the kind refuses what the block refuses and
 * nothing more. The four settings it bounds nowhere — the two hook-effect thresholds, the hook
 * read minimum, and the two Hill coefficients apart from their ordering — are checked only as
 * numbers here.
 */
const CONTRACT = {
  abundanceRef: check(isPlRef, "a reference to an abundance column"),
  concentrationColumnRef: check(isPlRef, "a reference to a concentration column"),
  binColumnRef: check(isPlRef, "a reference to a FACS bin column"),
  antigenColumnRef: check(isPlRef, "a reference to an antigen column"),
  sortFractionColumnRef: check(isPlRef, "a reference to a FACS sort fraction column"),
  targetAntigen: check(isString, "a string"),
  minReadsPerConcentration: check(isNumberFrom(1), "a number of 1 or more"),
  minConcentrationPoints: check(isNumberFrom(3), "a number of 3 or more"),
  r2ThresholdGood: check(isFraction, "a number between 0 and 1"),
  r2ThresholdFailed: check(isFraction, "a number between 0 and 1"),
  nMin: check(isFiniteNumber, "a number"),
  nMax: check(isFiniteNumber, "a number"),
  hookEffectThresholdBin: check(isFiniteNumber, "a number"),
  hookEffectThresholdNoBin: check(isFiniteNumber, "a number"),
  hookEffectMinReads: check(isFiniteNumber, "a number"),
  customBlockLabel: check(isString, "a string"),
} satisfies { [K in keyof Required<BlockParams>]: Check<NonNullable<BlockParams[K]>> };
