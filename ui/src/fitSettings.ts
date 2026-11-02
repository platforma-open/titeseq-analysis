import type { Equation } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { EQUATIONS } from "@platforma-open/platforma-open.titeseq-analysis.model";

// Exhaustive over the model's equation union: a new equation fails to compile until it is named.
const EQUATION_LABELS: Record<Equation, string> = {
  log_fluorescence: "Log fluorescence fit",
  hill_gate_index: "Mean bin fit",
};

export const EQUATION_OPTIONS: readonly { readonly value: Equation; readonly label: string }[] =
  EQUATIONS.map((value) => ({ value, label: EQUATION_LABELS[value] }));
