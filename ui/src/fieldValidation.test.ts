import { ref } from "vue";
import { describe, expect, it } from "vitest";
import { aboveZero, boundedBy, emptyError, wholeNumber } from "./fieldValidation";

describe("empty", () => {
  it("rejects an unset value", () => {
    expect(emptyError(undefined)).toBe("Value is required");
  });

  it.each([0, 0.4, 2 / 3, 20])("accepts %s", (v) => {
    expect(emptyError(v)).toBeUndefined();
  });
});

describe("whole number", () => {
  it.each([0, 1, 3, 20, -4])("accepts %s", (v) => {
    expect(wholeNumber(v)).toBeUndefined();
  });

  it.each([0.5, 3.5, 2.000001])("rejects %s", (v) => {
    expect(wholeNumber(v)).toBe("Value must be a whole number");
  });

  // Below the field's own `min-value`, which PlNumberField reports before `validate` runs.
  it("leaves the lower bound to the field", () => {
    expect(wholeNumber(-4)).toBeUndefined();
  });
});

describe("above zero", () => {
  // The one value `min-value="0"` admits that the plateau floor cannot use.
  it("rejects 0", () => {
    expect(aboveZero(0)).toBe("Value must be greater than 0");
  });

  it.each([Number.MIN_VALUE, 0.4, 1])("accepts %s", (v) => {
    expect(aboveZero(v)).toBeUndefined();
  });
});

describe("a bound that is another field's value", () => {
  const cases = [
    { relation: "≥", bound: 0.5, rejects: [0.4, 0], accepts: [0.5, 0.8] },
    { relation: "≤", bound: 0.8, rejects: [0.9, 1], accepts: [0.8, 0.5] },
    { relation: "<", bound: 2, rejects: [2, 3], accepts: [1.9, 0] },
    { relation: ">", bound: 0.5, rejects: [0.5, 0.1], accepts: [0.6, 2] },
  ] as const;

  it.each(cases)("$relation $bound", ({ relation, bound, rejects, accepts }) => {
    const check = boundedBy(relation, () => bound, "the other field");
    for (const v of rejects)
      expect(check(v), `${v}`).toBe(`Value must be ${relation} the other field`);
    for (const v of accepts) expect(check(v), `${v}`).toBeUndefined();
  });

  it.each(cases)("$relation passes while the other field is unset", ({ relation }) => {
    expect(boundedBy(relation, () => undefined, "the other field")(1)).toBeUndefined();
  });

  // PlNumberField calls `validate` inside a computed, so the verdict has to follow the other
  // field rather than the value it held when the drawer mounted.
  it("re-reads the other field on every call", () => {
    const nMax = ref<number | undefined>(2);
    const nMinBelowMax = boundedBy("<", () => nMax.value, "Hill coefficient — max");

    expect(nMinBelowMax(1.5)).toBeUndefined();
    nMax.value = 1;
    expect(nMinBelowMax(1.5)).toBe("Value must be < Hill coefficient — max");
  });
});

describe("the two drawer pairs disagree in exactly one direction", () => {
  // R² Good must be at or above Failed, so one message appears, on the field that is out of
  // place relative to the other — never both at once for a pair that is merely ordered wrong.
  it("R² thresholds", () => {
    const good = ref<number | undefined>(0.8);
    const failed = ref<number | undefined>(0.5);
    const goodCheck = boundedBy("≥", () => failed.value, "R² threshold (Failed)");
    const failedCheck = boundedBy("≤", () => good.value, "R² threshold (Good)");

    expect(goodCheck(good.value ?? 0)).toBeUndefined();
    expect(failedCheck(failed.value ?? 0)).toBeUndefined();

    good.value = 0.3;
    expect(goodCheck(0.3)).toBe("Value must be ≥ R² threshold (Failed)");
    expect(failedCheck(0.5)).toBe("Value must be ≤ R² threshold (Good)");
  });

  // nMin < nMax is strict: equal values leave no interval, and both fields say so.
  it("Hill coefficients", () => {
    const nMin = ref<number | undefined>(0.5);
    const nMax = ref<number | undefined>(2);
    const minCheck = boundedBy("<", () => nMax.value, "Hill coefficient — max");
    const maxCheck = boundedBy(">", () => nMin.value, "Hill coefficient — min");

    expect(minCheck(0.5)).toBeUndefined();
    expect(maxCheck(2)).toBeUndefined();

    nMin.value = 2;
    expect(minCheck(2)).toBe("Value must be < Hill coefficient — max");
    expect(maxCheck(2)).toBe("Value must be > Hill coefficient — min");
  });
});
