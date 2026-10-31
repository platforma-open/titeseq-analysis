import { describe, expect, it } from "vitest";
import { initBlockData } from "./index";
import type { FitSettings } from "./presets";
import { HILL_GATE_INDEX_DEFAULTS, PRESETS, matchesPreset, presetById } from "./presets";
import {
  SETTINGS_FILE_FORMAT,
  SETTINGS_FILE_VERSION,
  readSettingsFile,
  writeSettingsFile,
} from "./settingsFile";

const settings: FitSettings = {
  equation: "log_fluorescence",
  plateauFloor: 0.4,
  ...HILL_GATE_INDEX_DEFAULTS,
  nMax: 3,
};

/** A file body with every required field, so each test alters exactly one thing. */
const body = (patch: Record<string, unknown> = {}) =>
  JSON.stringify({
    format: SETTINGS_FILE_FORMAT,
    version: SETTINGS_FILE_VERSION,
    equation: "log_fluorescence",
    plateauFloor: 0.4,
    ...patch,
  });

describe("writing and reading a settings file", () => {
  it("is a round trip", () => {
    const read = readSettingsFile(writeSettingsFile(settings));
    expect(read).toEqual({ status: "ready", settings });
  });

  it("carries the fit settings and nothing that belongs to one project", () => {
    const written: unknown = JSON.parse(writeSettingsFile(settings));
    expect(Object.keys(written as object).sort()).toEqual(
      ["format", "version", ...Object.keys(settings)].sort(),
    );
  });

  it("drops column references a hand-edited file carries", () => {
    const read = readSettingsFile(
      body({ abundanceRef: { __isRef: true, blockId: "b1", name: "abundance" } }),
    );
    expect(read).toEqual({
      status: "ready",
      settings: { equation: "log_fluorescence", plateauFloor: 0.4 },
    });
  });
});

describe("refusing a file", () => {
  it("says so when the text is not JSON", () => {
    expect(readSettingsFile("not json").status).toBe("unparsable");
  });

  it.each([
    ["the format is another tool's", body({ format: "other/settings" }), "'format'"],
    ["the version is unknown", body({ version: 2 }), "'version'"],
    ["the equation is not one we fit", body({ equation: "poisson" }), "'equation'"],
    ["the plateau floor is out of range", body({ plateauFloor: 0 }), "'plateauFloor'"],
    ["a threshold is out of range", body({ r2ThresholdGood: 1.5 }), "'r2ThresholdGood'"],
    [
      "the equation is missing",
      JSON.stringify({ format: SETTINGS_FILE_FORMAT, version: 1 }),
      "'equation'",
    ],
  ])("refuses a file where %s, and says which field", (_name, text, field) => {
    const read = readSettingsFile(text);
    expect(read.status).toBe("rejected");
    if (read.status !== "rejected") throw new Error("unreachable");
    expect(read.reason).toContain(field);
  });

  it("takes the whole file or none of it", () => {
    // A good equation next to a bad threshold: neither is applied.
    expect(readSettingsFile(body({ r2ThresholdGood: 1.5 })).status).toBe("rejected");
  });
});

describe("the presets", () => {
  it("names each published analysis once", () => {
    expect(PRESETS.map((p) => p.id)).toEqual(["adams2016", "phillips2021", "starr"]);
  });

  it("writes the equation each analysis used", () => {
    expect(presetById("adams2016").settings.equation).toBe("log_fluorescence");
    expect(presetById("phillips2021").settings.plateauFloor).toBe(0.4);
    expect(presetById("starr").settings.equation).toBe("hill_gate_index");
  });

  it("reads back as applied once its settings are in place", () => {
    const preset = presetById("phillips2021");
    const applied: FitSettings = { ...initBlockData(), ...preset.settings };
    expect(matchesPreset(preset, applied)).toBe(true);
    expect(matchesPreset(preset, { ...applied, plateauFloor: 0.5 })).toBe(false);
  });

  it("holds the block's own defaults for the settings only the gate-index fit reads", () => {
    const fresh = initBlockData();
    for (const [field, value] of Object.entries(HILL_GATE_INDEX_DEFAULTS)) {
      expect(fresh[field as keyof typeof HILL_GATE_INDEX_DEFAULTS]).toBe(value);
    }
  });
});
