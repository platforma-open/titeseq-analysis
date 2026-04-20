import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';

import { uniquePlId } from '@platforma-sdk/model';
import type { ImportFileHandle } from '@platforma-sdk/model';
import type {
  BlockArgs as SamplesAndDataBlockArgs,
  PlId,
} from '@platforma-open/milaboratories.samples-and-data.model';

export type FixtureVariantName = 'bin_mode' | 'no_bin_mode' | 'antigen';

type RawVariant = {
  variant: FixtureVariantName;
  sample_ids: string[];
  concentration: Record<string, number>;
  bin?: Record<string, number>;
  antigen?: Record<string, string>;
  target_antigen?: string;
  distractor_clonotypes?: string[];
  clonotypes: Record<string, Record<string, unknown>>;
  cdr3_to_label: Record<string, string>;
};

export type FixtureManifest = {
  fixture_version: string;
  master_seed: number;
  concentrations_m: number[];
  bins: number[];
  mixcr_tsv_headers: string[];
  v_gene: string;
  j_gene: string;
  variants: Record<FixtureVariantName, RawVariant>;
};

export type LoadedVariant = {
  variant: FixtureVariantName;
  /** Stable fixture sample names (e.g. "c0_b1") → platforma sample ids. */
  sampleIdByName: Record<string, PlId>;
  /** PlId → label (same as fixture name). */
  sampleLabels: Record<PlId, string>;
  /** Metadata column ids keyed by purpose. */
  metadataIds: {
    concentration: PlId;
    bin?: PlId;
    antigen?: PlId;
  };
  /** Dataset id (S&D resource). */
  datasetId: PlId;
  /** The raw manifest entry for this variant (for test assertions). */
  raw: RawVariant;
};

const MANIFEST_PATH = resolve(__dirname, '..', 'fixtures', 'manifest.json');

export function loadManifest(): FixtureManifest {
  return JSON.parse(readFileSync(MANIFEST_PATH, 'utf8')) as FixtureManifest;
}

export function fixtureTsvPath(variant: FixtureVariantName, sampleName: string): string {
  return resolve(__dirname, '..', 'fixtures', 'data', variant, `sample_${sampleName}.tsv`);
}

export type Helpers = {
  getLocalFileHandle: (path: string) => Promise<string>;
};

/**
 * Build the S&D `BlockArgs` for a fixture variant and materialize local file
 * handles for every sample TSV. Returns both the args payload (to pass to
 * `project.setBlockArgs(sndBlockId, ...)`) and the PlId lookup tables the
 * tests need to cross-reference outputs back to fixture sample names.
 */
export async function prepareSamplesAndDataArgs(
  variant: FixtureVariantName,
  helpers: Helpers,
): Promise<{ args: SamplesAndDataBlockArgs; loaded: LoadedVariant }> {
  const manifest = loadManifest();
  const raw = manifest.variants[variant];
  if (!raw) throw new Error(`unknown fixture variant: ${variant}`);

  const sampleIdByName: Record<string, PlId> = {};
  const sampleLabels: Record<PlId, string> = {};
  for (const name of raw.sample_ids) {
    const id = uniquePlId() as unknown as PlId;
    sampleIdByName[name] = id;
    sampleLabels[id] = name;
  }

  // Import file handles for every sample TSV. S&D will copy these into the
  // project on run; the handles are only valid for this test session.
  const datasetData: Record<PlId, ImportFileHandle> = {};
  for (const name of raw.sample_ids) {
    datasetData[sampleIdByName[name]] = (await helpers.getLocalFileHandle(
      fixtureTsvPath(variant, name),
    )) as unknown as ImportFileHandle;
  }

  const concentrationId = uniquePlId() as unknown as PlId;
  const metadataIds: LoadedVariant['metadataIds'] = { concentration: concentrationId };

  const metadata: SamplesAndDataBlockArgs['metadata'] = [
    {
      id: concentrationId,
      label: 'antigen_conc_M',
      global: true,
      valueType: 'Double',
      data: Object.fromEntries(
        raw.sample_ids.map((n) => [sampleIdByName[n], raw.concentration[n]]),
      ),
    },
  ];

  if (raw.bin) {
    const binId = uniquePlId() as unknown as PlId;
    metadataIds.bin = binId;
    metadata.push({
      id: binId,
      label: 'bin_number',
      global: true,
      valueType: 'Long',
      data: Object.fromEntries(
        raw.sample_ids.map((n) => [sampleIdByName[n], raw.bin![n]]),
      ),
    });
  }

  if (raw.antigen) {
    const antigenId = uniquePlId() as unknown as PlId;
    metadataIds.antigen = antigenId;
    metadata.push({
      id: antigenId,
      label: 'antigen',
      global: true,
      valueType: 'String',
      data: Object.fromEntries(
        raw.sample_ids.map((n) => [sampleIdByName[n], raw.antigen![n]]),
      ),
    });
  }

  const datasetId = uniquePlId() as unknown as PlId;

  const args: SamplesAndDataBlockArgs = {
    metadata,
    sampleIds: raw.sample_ids.map((n) => sampleIdByName[n]),
    sampleLabelColumnLabel: 'Sample Name',
    sampleLabels,
    datasets: [
      {
        id: datasetId,
        label: `Titeseq ${variant} fixture`,
        content: {
          type: 'Xsv',
          xsvType: 'tsv',
          gzipped: false,
          data: datasetData,
        },
      },
    ],
  };

  return {
    args,
    loaded: { variant, sampleIdByName, sampleLabels, metadataIds, datasetId, raw },
  };
}
