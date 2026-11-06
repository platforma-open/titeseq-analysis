import type { AxisSpec, PColumnSpec, PObjectSpec } from "@platforma-sdk/model";
import { isPColumnSpec } from "@platforma-sdk/model";

export const AA_TO_NT = "pl7.app/repertoire/aaToNt";
export const MUTATION_COUNT = "pl7.app/repertoire/mutationCount";
export const PARENT_LINK = "pl7.app/repertoire/parentLink";
const PARENT_ID_AXIS = "pl7.app/repertoire/parentId";
const VARIANT_KEY_AXIS = "pl7.app/variantKey";
const ALPHABET = "pl7.app/alphabet";

/** What one row of the abundance is. Labels read it; it changes no computation. */
export type InputKind = "vdj-clonotype" | "variant" | "peptide";

export type ParentInput = typeof MUTATION_COUNT | typeof PARENT_LINK;

type UnidentifiableParent = {
  parentId: string;
  reason: "no-unmutated-variant" | "several-unmutated-variants";
};

/** The run summary the software writes in synonymous mode. */
export type SynonymousSummary = {
  ntVariants: number;
  aaVariants: number;
  unmappedNtVariants: number;
  ambiguousNtVariants: number;
  parentMode: boolean;
  parents?: number;
  unidentifiableParents?: UnidentifiableParent[];
  multiParentVariants?: number;
  /** Present when the nucleotide mutation count resolved. */
  ntParentMode?: boolean;
  ntParents?: number;
  ntUnidentifiableParents?: UnidentifiableParent[];
  ntMultiParentVariants?: number;
  /** Present when the amino-acid productive column resolved. */
  stopAaVariants?: number;
  /** The SD of log10 Kd,app over measured synonymous copies, per amino-acid variant with two or more. */
  synonymSpread?: { variants: number; median: number | null };
};

/**
 * `variantKey` is shared by the profiler, peptide extraction and imported antibody sets, so the
 * axis domain decides, not the name. A declared modality wins over the run-id keys.
 */
export function inputKindOf(axis: AxisSpec): InputKind {
  if (axis.name !== VARIANT_KEY_AXIS) return "vdj-clonotype";
  const domain = axis.domain ?? {};
  if (domain["pl7.app/modality"] === "vdj" || domain["pl7.app/vdj/clonotypingRunId"] !== undefined)
    return "vdj-clonotype";
  if (domain["pl7.app/peptide/extractionRunId"] !== undefined) return "peptide";
  return "variant";
}

/** The noun labels use for one row of the input. */
export function nounOf(kind: InputKind | undefined): "clonotype" | "variant" | "peptide" {
  switch (kind) {
    case "variant":
      return "variant";
    case "peptide":
      return "peptide";
    case "vdj-clonotype":
    case undefined:
      return "clonotype";
    default:
      throw new Error(`Unknown input kind: ${kind satisfies never}`);
  }
}

const domainKey = (axis: AxisSpec) =>
  JSON.stringify(Object.entries(axis.domain ?? {}).sort(([a], [b]) => a.localeCompare(b)));

const sameAxis = (a: AxisSpec, b: AxisSpec) =>
  a.name === b.name && a.type === b.type && domainKey(a) === domainKey(b);

const isLinker = (spec: PColumnSpec) => spec.annotations?.["pl7.app/isLinkerColumn"] === "true";

/** The linker's amino-acid axis when `spec` is an nt -> aa linker on `ntAxis`. */
const aaAxisOfLinker = (spec: PColumnSpec, ntAxis: AxisSpec): AxisSpec | undefined => {
  if (spec.name !== AA_TO_NT || !isLinker(spec) || spec.axesSpec.length !== 2) return undefined;
  const [a, b] = spec.axesSpec;
  const other = sameAxis(a, ntAxis) ? b : sameAxis(b, ntAxis) ? a : undefined;
  if (other?.name !== ntAxis.name || other.domain?.[ALPHABET] !== "aminoacid") return undefined;
  return other;
};

const isMutationCountOn = (spec: PColumnSpec, axis: AxisSpec) =>
  spec.name === MUTATION_COUNT && spec.axesSpec.length === 1 && sameAxis(spec.axesSpec[0], axis);

const isParentLinkOn = (spec: PColumnSpec, axis: AxisSpec) =>
  spec.name === PARENT_LINK &&
  isLinker(spec) &&
  spec.axesSpec.length === 2 &&
  spec.axesSpec.some((a) => sameAxis(a, axis)) &&
  spec.axesSpec.some((a) => a.name === PARENT_ID_AXIS);

/**
 * Mirrors the workflow's discovery: synonymous mode needs a nucleotide axis and exactly one linker;
 * several matches count as none. Returns the linker's amino-acid axis, or undefined when the mode is off.
 */
export function synonymousAaAxis(specs: PObjectSpec[], ntAxis: AxisSpec): AxisSpec | undefined {
  if (ntAxis.domain?.[ALPHABET] !== "nucleotide") return undefined;
  const found = specs
    .filter(isPColumnSpec)
    .map((s) => aaAxisOfLinker(s, ntAxis))
    .filter((a): a is AxisSpec => a !== undefined);
  return found.length === 1 ? found[0] : undefined;
}

/** The one amino-acid parent column found, when the other is missing; the delta is then off. */
export function halfMatchedParentInput(
  specs: PObjectSpec[],
  aaAxis: AxisSpec,
): ParentInput | undefined {
  const columns = specs.filter(isPColumnSpec);
  const mutationCount = columns.filter((s) => isMutationCountOn(s, aaAxis)).length === 1;
  const parentLink = columns.filter((s) => isParentLinkOn(s, aaAxis)).length === 1;
  if (mutationCount === parentLink) return undefined;
  return mutationCount ? MUTATION_COUNT : PARENT_LINK;
}

const PARENT_INPUT_LABELS: Record<ParentInput, string> = {
  [MUTATION_COUNT]: "amino-acid mutation count",
  [PARENT_LINK]: "amino-acid parent link",
};

const UNIDENTIFIABLE_REASONS: Record<UnidentifiableParent["reason"], string> = {
  "no-unmutated-variant": "no amino-acid variant has zero mutations",
  "several-unmutated-variants": "several amino-acid variants have zero mutations",
};

const NT_UNIDENTIFIABLE_REASONS: Record<UnidentifiableParent["reason"], string> = {
  "no-unmutated-variant": "no nucleotide variant has zero mutations",
  "several-unmutated-variants": "several nucleotide variants have zero mutations",
};

/** The messages the results pages show about the synonymous combination. */
export function synonymousMessages(
  halfMatched: ParentInput | undefined,
  summary: SynonymousSummary | undefined,
): string[] {
  const messages: string[] = [];
  if (halfMatched !== undefined) {
    const missing = halfMatched === MUTATION_COUNT ? PARENT_LINK : MUTATION_COUNT;
    messages.push(
      `The ${PARENT_INPUT_LABELS[halfMatched]} column was found, but the ` +
        `${PARENT_INPUT_LABELS[missing]} column was not. Kd,app is not compared with the parent.`,
    );
  }
  if (summary === undefined) return messages;
  const excluded = summary.unmappedNtVariants + summary.ambiguousNtVariants;
  if (excluded > 0)
    messages.push(
      `${excluded} of ${summary.ntVariants} nucleotide variants are left out of the amino-acid ` +
        `values: ${summary.unmappedNtVariants} translate to no amino-acid variant, ` +
        `${summary.ambiguousNtVariants} to several.`,
    );
  for (const p of summary.unidentifiableParents ?? [])
    messages.push(
      `Parent ${p.parentId} cannot be identified: ${UNIDENTIFIABLE_REASONS[p.reason]}. ` +
        "Its variants have no parent delta.",
    );
  if ((summary.multiParentVariants ?? 0) > 0)
    messages.push(
      `${summary.multiParentVariants} amino-acid variants link to several parents and have no parent delta.`,
    );
  for (const p of summary.ntUnidentifiableParents ?? [])
    messages.push(
      `Parent ${p.parentId} cannot be identified at nucleotide level: ` +
        `${NT_UNIDENTIFIABLE_REASONS[p.reason]}. Its nucleotide variants have no parent delta.`,
    );
  if ((summary.ntMultiParentVariants ?? 0) > 0)
    messages.push(
      `${summary.ntMultiParentVariants} nucleotide variants link to several parents and have no parent delta.`,
    );
  if ((summary.stopAaVariants ?? 0) > 0)
    messages.push(
      `${summary.stopAaVariants} amino-acid variants carry a stop codon. They have no Kd,app and are ` +
        "rated Failed; their nucleotide variants keep their own fits.",
    );
  return messages;
}
