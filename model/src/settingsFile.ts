import type { BlockParams } from "@platforma-open/platforma-open.titeseq-analysis.kind";
import { kind } from "@platforma-open/platforma-open.titeseq-analysis.kind";
import type { FitSettings } from "./presets";

export const SETTINGS_FILE_FORMAT = "titeseq-analysis/settings";
export const SETTINGS_FILE_VERSION = 1;
export const SETTINGS_FILE_NAME = "titeseq-settings.json";

/**
 * What a settings file carries. The equation and its plateau floor are always written and always
 * required; the settings only `hill_gate_index` reads are written too, but a file that omits one
 * leaves the block's own value standing.
 */
export type ImportedFitSettings = Pick<FitSettings, "equation" | "plateauFloor"> &
  Partial<Omit<FitSettings, "equation" | "plateauFloor">>;

/**
 * Every state a picked settings file can be in, as one value the UI renders without a cast.
 *
 * A file is read and checked here, in the model: the derivation is pure, the runtime memoizes it,
 * and the UI needs no watcher on a server-derived value — it shows this and writes to `data` only
 * when the scientist presses Apply.
 */
export type SettingsImport =
  /** No file picked. */
  | { readonly status: "none" }
  /** A file is picked; its bytes have not arrived yet. */
  | { readonly status: "loading" }
  /** The file is not JSON. */
  | { readonly status: "unparsable"; readonly reason: string }
  /** The file is JSON, but not a settings file this block applies. */
  | { readonly status: "rejected"; readonly reason: string }
  /** Every field present was read and is in range. */
  | { readonly status: "ready"; readonly settings: ImportedFitSettings };

const isRecord = (v: unknown): v is Record<string, unknown> =>
  typeof v === "object" && v !== null && !Array.isArray(v);

/** The settings, and only the settings, out of a block's state. */
export function writeSettingsFile(settings: FitSettings): string {
  return `${JSON.stringify(
    {
      format: SETTINGS_FILE_FORMAT,
      version: SETTINGS_FILE_VERSION,
      equation: settings.equation,
      plateauFloor: settings.plateauFloor,
      minReadsPerConcentration: settings.minReadsPerConcentration,
      minConcentrationPoints: settings.minConcentrationPoints,
      r2ThresholdGood: settings.r2ThresholdGood,
      r2ThresholdFailed: settings.r2ThresholdFailed,
      nMin: settings.nMin,
      nMax: settings.nMax,
      hookEffectThresholdBin: settings.hookEffectThresholdBin,
      hookEffectThresholdNoBin: settings.hookEffectThresholdNoBin,
      hookEffectMinReads: settings.hookEffectMinReads,
    },
    null,
    2,
  )}\n`;
}

/**
 * Read a settings file into the domain type, or say precisely why it was refused. The whole file
 * is taken or the whole file is refused: the first failure ends the read.
 *
 * Field names, types and bounds are the kind's own params contract, so a file and a project
 * template are checked by one piece of code. What the file may carry is narrower than that
 * contract: column references and labels belong to the project that wrote them, so they are
 * dropped rather than refused.
 */
export function readSettingsFile(text: string): SettingsImport {
  let body: unknown;
  try {
    body = JSON.parse(text);
  } catch (e) {
    return { status: "unparsable", reason: e instanceof Error ? e.message : String(e) };
  }

  if (!isRecord(body)) {
    return { status: "rejected", reason: "The file does not hold a settings object." };
  }
  const record = body;

  if (record["format"] !== SETTINGS_FILE_FORMAT) {
    return {
      status: "rejected",
      reason: `'format' must be "${SETTINGS_FILE_FORMAT}", not ${JSON.stringify(record["format"])}.`,
    };
  }
  if (record["version"] !== SETTINGS_FILE_VERSION) {
    return {
      status: "rejected",
      reason: `'version' must be ${SETTINGS_FILE_VERSION}, not ${JSON.stringify(record["version"])}.`,
    };
  }

  let params: BlockParams;
  try {
    params = kind.parseInitializationParams(record);
  } catch (e) {
    return { status: "rejected", reason: e instanceof Error ? e.message : String(e) };
  }

  const { equation, plateauFloor } = params;
  if (equation === undefined) return { status: "rejected", reason: "'equation' is required." };
  if (plateauFloor === undefined) {
    return { status: "rejected", reason: "'plateauFloor' is required." };
  }

  const settings: ImportedFitSettings = { equation, plateauFloor };
  const optional = [
    "minReadsPerConcentration",
    "minConcentrationPoints",
    "r2ThresholdGood",
    "r2ThresholdFailed",
    "nMin",
    "nMax",
    "hookEffectThresholdBin",
    "hookEffectThresholdNoBin",
    "hookEffectMinReads",
  ] as const;
  for (const field of optional) {
    const value = params[field];
    if (value !== undefined) settings[field] = value;
  }

  return { status: "ready", settings };
}
