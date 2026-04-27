import { computed } from "vue";
import { useApp } from "../app";

// Per-field error messages bound to each PlNumberField's `:error-message`
// prop. Inline only — page-level alerts are reserved for spec-based checks
// (concentration column label, sort-fraction without bin) emitted by the
// model's validationWarnings output and rendered in TiteseqPage.
//
// Each field returns one consolidated message covering the full valid range
// (e.g. "Must be between 0 and 1, and ≥ R² threshold (Failed)") instead of
// separate alerts per bound. For typed non-numeric input the message stays
// undefined so PlNumberField's own "Value is not a number" parse error wins.

const REQUIRED = "Value is required";

const isMissing = (v: number | undefined | null): boolean =>
  v === undefined || v === null || Number.isNaN(v);

// Integer ≥ min. One message covers non-integer + below-min + (when needed)
// the user-facing description "whole number ≥ N".
const validateInteger = (v: number | undefined | null, min: number): string | undefined => {
  if (isMissing(v)) return REQUIRED;
  if (!Number.isInteger(v as number) || (v as number) < min) {
    return `Must be a whole number ≥ ${min}`;
  }
  return undefined;
};

// Float in [0, 1] with an optional cross-field bound. `relation` is "≥" when
// the bound is a lower limit (this field must be ≥ bound) and "≤" when the
// bound is an upper limit. The bound is ignored if not a finite number.
const validate01WithBound = (
  v: number | undefined | null,
  bound: number | undefined | null,
  relation: "≥" | "≤",
  boundLabel: string,
): string | undefined => {
  if (isMissing(v)) return REQUIRED;
  const inRange = (v as number) >= 0 && (v as number) <= 1;
  const passesBound =
    isMissing(bound) ||
    (relation === "≥" ? (v as number) >= (bound as number) : (v as number) <= (bound as number));
  if (!inRange || !passesBound) {
    return `Must be between 0 and 1, and ${relation} ${boundLabel}`;
  }
  return undefined;
};

// Float ≥ 0 with an optional strict cross-field bound. Used for nMin/nMax —
// the spec requires nMin < nMax (strict) to keep the in-range interval
// non-empty.
const validateNonNegWithStrict = (
  v: number | undefined | null,
  bound: number | undefined | null,
  relation: "<" | ">",
  boundLabel: string,
): string | undefined => {
  if (isMissing(v)) return REQUIRED;
  const nonNeg = (v as number) >= 0;
  const passesBound =
    isMissing(bound) ||
    (relation === "<" ? (v as number) < (bound as number) : (v as number) > (bound as number));
  if (!nonNeg || !passesBound) {
    return `Must be ≥ 0 and ${relation} ${boundLabel}`;
  }
  return undefined;
};

// Plain non-negative float, no cross-field bound (hookEffectThreshold*).
const validateNonNeg = (v: number | undefined | null): string | undefined => {
  if (isMissing(v)) return REQUIRED;
  if ((v as number) < 0) return "Must be ≥ 0";
  return undefined;
};

export function useFieldValidation() {
  const app = useApp();

  // Per-field reactive errors. Each computed produces at most one message —
  // either "Value is required" (when undefined/NaN) or a single consolidated
  // bounds message that covers every range/cross-field violation.
  const minReadsError = computed(() => validateInteger(app.model.data.minReadsPerConcentration, 1));
  const minConcPointsError = computed(() =>
    validateInteger(app.model.data.minConcentrationPoints, 3),
  );
  const r2GoodError = computed(() =>
    validate01WithBound(
      app.model.data.r2ThresholdGood,
      app.model.data.r2ThresholdFailed,
      "≥",
      "R² threshold (Failed)",
    ),
  );
  const r2FailedError = computed(() =>
    validate01WithBound(
      app.model.data.r2ThresholdFailed,
      app.model.data.r2ThresholdGood,
      "≤",
      "R² threshold (Good)",
    ),
  );
  const nMinError = computed(() =>
    validateNonNegWithStrict(
      app.model.data.nMin,
      app.model.data.nMax,
      "<",
      "Hill coefficient — max",
    ),
  );
  const nMaxError = computed(() =>
    validateNonNegWithStrict(
      app.model.data.nMax,
      app.model.data.nMin,
      ">",
      "Hill coefficient — min",
    ),
  );
  const hookThresholdBinError = computed(() =>
    validateNonNeg(app.model.data.hookEffectThresholdBin),
  );
  const hookThresholdNoBinError = computed(() =>
    validateNonNeg(app.model.data.hookEffectThresholdNoBin),
  );
  const hookMinReadsError = computed(() => validateInteger(app.model.data.hookEffectMinReads, 0));

  return {
    minReadsError,
    minConcPointsError,
    r2GoodError,
    r2FailedError,
    nMinError,
    nMaxError,
    hookThresholdBinError,
    hookThresholdNoBinError,
    hookMinReadsError,
  };
}
