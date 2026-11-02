import type {
  Equation,
  ReplicateErrorModel,
} from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  EQUATIONS,
  REPLICATE_ERROR_MODELS,
} from "@platforma-open/platforma-open.titeseq-analysis.model";

// Exhaustive over the model's equation union: a new equation fails to compile until it is named.
const EQUATION_LABELS: Record<Equation, string> = {
  log_fluorescence: "Log fluorescence fit",
  hill_gate_index: "Mean bin fit",
};

export const EQUATION_OPTIONS: readonly { readonly value: Equation; readonly label: string }[] =
  EQUATIONS.map((value) => ({ value, label: EQUATION_LABELS[value] }));

const REPLICATE_ERROR_MODEL_LABELS: Record<ReplicateErrorModel, string> = {
  sd_over_sqrt_n: "SD / √n",
  phillips: "Phillips 2021",
};

export const REPLICATE_ERROR_MODEL_OPTIONS: readonly {
  readonly value: ReplicateErrorModel;
  readonly label: string;
}[] = REPLICATE_ERROR_MODELS.map((value) => ({
  value,
  label: REPLICATE_ERROR_MODEL_LABELS[value],
}));
