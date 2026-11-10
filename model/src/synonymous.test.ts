import type { AxisSpec, PColumnSpec } from "@platforma-sdk/model";
import { describe, expect, it } from "vitest";
import type { SynonymousSummary } from "./synonymous";
import {
  MUTATION_COUNT,
  PARENT_LINK,
  halfMatchedParentInput,
  inputKindOf,
  nounOf,
  synonymousAaAxis,
  synonymousMessages,
} from "./synonymous";

const variantAxis = (alphabet: string, run = "run1"): AxisSpec => ({
  name: "pl7.app/variantKey",
  type: "String",
  domain: {
    "pl7.app/repertoire/extractionRunId": run,
    "pl7.app/modality": "amplicon",
    "pl7.app/alphabet": alphabet,
  },
  annotations: { "pl7.app/label": alphabet },
});
const parentAxis: AxisSpec = { name: "pl7.app/repertoire/parentId", type: "String" };

const column = (name: string, axes: AxisSpec[], linker = false): PColumnSpec => ({
  kind: "PColumn",
  name,
  valueType: "Int",
  axesSpec: axes,
  annotations: linker ? { "pl7.app/isLinkerColumn": "true" } : {},
});

const nt = variantAxis("nucleotide");
const aa = variantAxis("aminoacid");
const linker = column("pl7.app/repertoire/aaToNt", [aa, nt], true);

describe("inputKindOf", () => {
  it.for([
    [{ name: "pl7.app/vdj/clonotypeKey", type: "String" }, "vdj-clonotype"],
    [{ name: "pl7.app/vdj/scClonotypeKey", type: "String" }, "vdj-clonotype"],
    [nt, "variant"],
    [aa, "variant"],
    [
      { name: "pl7.app/variantKey", type: "String", domain: { "pl7.app/modality": "vdj" } },
      "vdj-clonotype",
    ],
    [
      {
        name: "pl7.app/variantKey",
        type: "String",
        domain: { "pl7.app/vdj/clonotypingRunId": "r" },
      },
      "vdj-clonotype",
    ],
    [
      {
        name: "pl7.app/variantKey",
        type: "String",
        domain: { "pl7.app/peptide/extractionRunId": "p" },
      },
      "peptide",
    ],
  ] as [AxisSpec, string][])("%o is %s", ([axis, kind], { expect }) => {
    expect(inputKindOf(axis)).toBe(kind);
  });

  it("names a row by kind, and a clonotype before any input", () => {
    expect(nounOf("variant")).toBe("variant");
    expect(nounOf("peptide")).toBe("peptide");
    expect(nounOf("vdj-clonotype")).toBe("clonotype");
    expect(nounOf(undefined)).toBe("clonotype");
  });
});

describe("synonymousAaAxis", () => {
  it("is the linker's other axis on a nucleotide dataset", () => {
    expect(synonymousAaAxis([column("x", [nt]), linker], nt)).toEqual(aa);
  });

  it("is off on an amino-acid dataset", () => {
    expect(synonymousAaAxis([linker], aa)).toBeUndefined();
  });

  it("is off with no linker or with several", () => {
    expect(synonymousAaAxis([], nt)).toBeUndefined();
    const other = column("pl7.app/repertoire/aaToNt", [variantAxis("aminoacid", "run2"), nt], true);
    expect(synonymousAaAxis([linker, other], nt)).toBeUndefined();
  });

  it("needs the linker annotation", () => {
    expect(synonymousAaAxis([column("pl7.app/repertoire/aaToNt", [aa, nt])], nt)).toBeUndefined();
  });
});

describe("halfMatchedParentInput", () => {
  // The workflow's query requires the column's own alphabet domain, so the model does too.
  const mutationCount = {
    ...column(MUTATION_COUNT, [aa]),
    domain: { "pl7.app/alphabet": "aminoacid" },
  };
  const parentLink = column(PARENT_LINK, [parentAxis, aa], true);

  it("is silent when both or neither resolve", () => {
    expect(halfMatchedParentInput([mutationCount, parentLink], aa)).toBeUndefined();
    expect(halfMatchedParentInput([], aa)).toBeUndefined();
  });

  it("does not count a mutation count without the amino-acid column domain", () => {
    expect(halfMatchedParentInput([column(MUTATION_COUNT, [aa]), parentLink], aa)).toBe(
      PARENT_LINK,
    );
  });

  it("names the column found when the other is missing", () => {
    expect(halfMatchedParentInput([mutationCount], aa)).toBe(MUTATION_COUNT);
    expect(halfMatchedParentInput([parentLink], aa)).toBe(PARENT_LINK);
  });

  it("reads the nucleotide grain's parent columns as absent", () => {
    expect(halfMatchedParentInput([column(MUTATION_COUNT, [nt]), parentLink], aa)).toBe(
      PARENT_LINK,
    );
  });
});

describe("synonymousMessages", () => {
  const summary: SynonymousSummary = {
    ntVariants: 10,
    aaVariants: 4,
    unmappedNtVariants: 0,
    ambiguousNtVariants: 0,
    parentMode: true,
    parents: 1,
    unidentifiableParents: [],
    multiParentVariants: 0,
  };

  it("says nothing for a clean run", () => {
    expect(synonymousMessages(undefined, summary)).toEqual([]);
    expect(synonymousMessages(undefined, undefined)).toEqual([]);
  });

  it("names the missing parent column", () => {
    expect(synonymousMessages(MUTATION_COUNT, undefined)).toEqual([
      "The amino-acid mutation count column was found, but the amino-acid parent link column " +
        "was not. Kd,app is not compared with the parent.",
    ]);
  });

  it("counts the nucleotide variants left out", () => {
    const [message] = synonymousMessages(undefined, {
      ...summary,
      unmappedNtVariants: 2,
      ambiguousNtVariants: 1,
    });
    expect(message).toBe(
      "3 of 10 nucleotide variants are left out of the amino-acid values: " +
        "2 translate to no amino-acid variant, 1 to several.",
    );
  });

  it("names an nt parent that cannot be identified", () => {
    const messages = synonymousMessages(undefined, {
      ...summary,
      ntParentMode: true,
      ntParents: 2,
      ntUnidentifiableParents: [{ parentId: "P2", reason: "several-unmutated-variants" }],
      ntMultiParentVariants: 4,
    });
    expect(messages).toEqual([
      "Parent P2 cannot be identified at nucleotide level: several nucleotide variants have zero " +
        "mutations. Its nucleotide variants have no parent delta.",
      "4 nucleotide variants link to several parents and have no parent delta.",
    ]);
  });

  it("counts the stop variants", () => {
    expect(synonymousMessages(undefined, { ...summary, stopAaVariants: 3 })).toEqual([
      "3 amino-acid variants carry a stop codon. They have no Kd,app and are rated Failed; " +
        "their nucleotide variants keep their own fits.",
    ]);
    expect(synonymousMessages(undefined, { ...summary, stopAaVariants: 0 })).toEqual([]);
  });

  it("names a parent that cannot be identified", () => {
    const messages = synonymousMessages(undefined, {
      ...summary,
      unidentifiableParents: [{ parentId: "P1", reason: "no-unmutated-variant" }],
      multiParentVariants: 2,
    });
    expect(messages).toEqual([
      "Parent P1 cannot be identified: no amino-acid variant has zero mutations. " +
        "Its variants have no parent delta.",
      "2 amino-acid variants link to several parents and have no parent delta.",
    ]);
  });
});
