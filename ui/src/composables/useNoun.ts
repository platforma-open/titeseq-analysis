import { nounOf } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { computed } from "vue";
import { useApp } from "../app";

/** The word user-visible text uses for one row of the input: clonotype, variant or peptide. */
export function useNoun() {
  const app = useApp();
  const noun = computed(() => nounOf(app.model.outputs.inputKind));
  const Noun = computed(() => noun.value.charAt(0).toUpperCase() + noun.value.slice(1));
  return { noun, Noun };
}
