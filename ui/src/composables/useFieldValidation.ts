import type { ValidationIssue } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

// Single source of truth for value-based field validation. Lives in the UI
// layer so error messages update synchronously as the user types — the model's
// validationWarnings output also emits these but lags behind the
// mutation→server→patch round-trip, which made the page-level alert appear
// "stuck" after entering a valid value.
//
// The model's validationWarnings is still authoritative for spec-based checks
// (concentration column label, antigen/target pairing) that need ctx.resultPool;
// TiteseqPage merges both so the user sees everything.

const validateIntField = (v: number | undefined | null, min: number): string | undefined => {
  if (v === undefined || v === null || Number.isNaN(v)) return "Value is required";
  if (!Number.isInteger(v)) return "Must be a whole number";
  if (v < min) return `Must be ≥ ${min}`;
  return undefined;
};

const validateFloatField = (
  v: number | undefined | null,
  opts: { min?: number; max?: number } = {},
): string | undefined => {
  if (v === undefined || v === null || Number.isNaN(v)) return "Value is required";
  if (opts.min !== undefined && v < opts.min) return `Must be ≥ ${opts.min}`;
  if (opts.max !== undefined && v > opts.max) return `Must be ≤ ${opts.max}`;
  return undefined;
};

export function useFieldValidation() {
  const app = useApp();

  const isBinMode = computed(() => app.model.data.binColumnRef !== undefined);

  // Per-field reactive errors. Pass these to PlNumberField via :error-message
  // to force the red contour and inline error text on the field itself.
  const minReadsError = computed(() =>
    validateIntField(app.model.data.minReadsPerConcentration, 1),
  );
  const minConcPointsError = computed(() =>
    validateIntField(app.model.data.minConcentrationPoints, 3),
  );
  const r2GoodError = computed(() =>
    validateFloatField(app.model.data.r2ThresholdGood, { min: 0, max: 1 }),
  );
  const r2FailedError = computed(() =>
    validateFloatField(app.model.data.r2ThresholdFailed, { min: 0, max: 1 }),
  );
  const nMinError = computed(() => validateFloatField(app.model.data.nMin, { min: 0 }));
  const nMaxError = computed(() => validateFloatField(app.model.data.nMax, { min: 0 }));
  const hookThresholdBinError = computed(() =>
    validateFloatField(app.model.data.hookEffectThresholdBin, { min: 0 }),
  );
  const hookThresholdNoBinError = computed(() =>
    validateFloatField(app.model.data.hookEffectThresholdNoBin, { min: 0 }),
  );
  const hookMinReadsError = computed(() => validateIntField(app.model.data.hookEffectMinReads, 0));

  // Aggregated page-level alerts. Includes per-field errors with their human
  // label, plus cross-field invariant checks (r2Failed > r2Good; nMin >= nMax).
  // Conditional fields (hookEffectThreshold*) only contribute the entry for
  // the active mode — the other field is hidden in the UI.
  const warnings = computed<ValidationIssue[]>(() => {
    const issues: ValidationIssue[] = [];
    const data = app.model.data;

    const fieldChecks: Array<[string, string | undefined]> = [
      ["Min reads per concentration", minReadsError.value],
      ["Min concentration points", minConcPointsError.value],
      ["R² threshold (Good)", r2GoodError.value],
      ["R² threshold (Failed)", r2FailedError.value],
      ["Hill coefficient nMin", nMinError.value],
      ["Hill coefficient nMax", nMaxError.value],
      isBinMode.value
        ? ["Hook effect signal-drop threshold (bin mode)", hookThresholdBinError.value]
        : ["Hook effect signal-drop threshold (frequency mode)", hookThresholdNoBinError.value],
      ["Min reads for hook check", hookMinReadsError.value],
    ];

    for (const [label, err] of fieldChecks) {
      if (err) {
        // "Value is required" → "<label> is required."; otherwise "<label>: <message>."
        const message = err === "Value is required" ? `${label} is required.` : `${label}: ${err}.`;
        issues.push({ severity: "error", message });
      }
    }

    if (
      data.r2ThresholdFailed !== undefined &&
      data.r2ThresholdGood !== undefined &&
      data.r2ThresholdFailed > data.r2ThresholdGood
    ) {
      issues.push({
        severity: "error",
        message: "Failed R² threshold must be ≤ Good R² threshold.",
      });
    }

    if (data.nMin !== undefined && data.nMax !== undefined && data.nMin >= data.nMax) {
      issues.push({
        severity: "error",
        message: "Hill coefficient nMin must be strictly less than nMax.",
      });
    }

    return issues;
  });

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
    warnings,
  };
}
