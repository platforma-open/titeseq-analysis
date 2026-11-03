// The checks the settings drawer's number fields take.
//
// PlNumberField's `min-value`/`max-value` cover the plain bounds with the wording every block
// shows. A `validate` callback here covers only what those props cannot express — a whole number,
// an exclusive bound, and the two pairs whose bound is the other field's value.
//
// PlNumberField calls `validate` inside a computed, so a callback that reads another setting
// re-runs when that setting changes.

/**
 * The fields are not `required`: a required PlNumberField keeps its model value when emptied and
 * its clear button must restore a number. Cleared, the model is unset; this marks the field red
 * and `args()` disables Run.
 */
export const emptyError = (v: number | undefined): string | undefined =>
  v === undefined ? "Value is required" : undefined;

/** Counts and read floors. `min-value` admits 3.5 where the fit needs 3. */
export const wholeNumber = (v: number): string | undefined =>
  Number.isInteger(v) ? undefined : "Value must be a whole number";

/** `min-value="0"` admits 0, and a plateau floor of 0 asks the curve to rise by nothing. */
export const aboveZero = (v: number): string | undefined =>
  v > 0 ? undefined : "Value must be greater than 0";

type Relation = "≥" | "≤" | "<" | ">";

const FAILS: Record<Relation, (v: number, bound: number) => boolean> = {
  "≥": (v, bound) => v < bound,
  "≤": (v, bound) => v > bound,
  "<": (v, bound) => v >= bound,
  ">": (v, bound) => v <= bound,
};

/**
 * `v` must stand in `relation` to another field's value, read through `bound` on every call. An
 * unset counterpart passes: a pair is half-set for as long as it takes to type the second value.
 */
export const boundedBy =
  (relation: Relation, bound: () => number | undefined, label: string) =>
  (v: number): string | undefined => {
    const b = bound();
    return b !== undefined && FAILS[relation](v, b)
      ? `Value must be ${relation} ${label}`
      : undefined;
  };
