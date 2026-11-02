import type {
  Equation,
  ReplicateErrorModel,
  SynonymousWeighting,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  EQUATIONS,
  REPLICATE_ERROR_MODELS,
  SYNONYMOUS_WEIGHTINGS,
} from "@platforma-open/platforma-open.titeseq-analysis.model";

// Exhaustive over the model's equation union: a new equation fails to compile until it is named.
const EQUATION_LABELS: Record<Equation, string> = {
  log_fluorescence: "Log fluorescence fit",
  hill_gate_index: "Mean bin fit",
};

export const EQUATION_OPTIONS: readonly { readonly value: Equation; readonly label: string }[] =
  EQUATIONS.map((value) => ({ value, label: EQUATION_LABELS[value] }));

// Named by their formula: the two differ only in which SD and which denominator.
const REPLICATE_ERROR_MODEL_LABELS: Record<ReplicateErrorModel, string> = {
  sd_over_sqrt_n: "Sample SD / √n",
  pop_sd_over_sqrt_n_plus_1: "Population SD / √(n+1)",
};

export const REPLICATE_ERROR_MODEL_OPTIONS: readonly {
  readonly value: ReplicateErrorModel;
  readonly label: string;
}[] = REPLICATE_ERROR_MODELS.map((value) => ({
  value,
  label: REPLICATE_ERROR_MODEL_LABELS[value],
}));

const SYNONYMOUS_WEIGHTING_LABELS: Record<SynonymousWeighting, string> = {
  none: "Plain mean",
  reads: "Weighted by reads",
};

export const SYNONYMOUS_WEIGHTING_OPTIONS: readonly {
  readonly value: SynonymousWeighting;
  readonly label: string;
}[] = SYNONYMOUS_WEIGHTINGS.map((value) => ({ value, label: SYNONYMOUS_WEIGHTING_LABELS[value] }));
