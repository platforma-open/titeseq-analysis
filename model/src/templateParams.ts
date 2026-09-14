import type { BlockParams } from "@platforma-open/platforma-open.titeseq-analysis.kind";
import type { BlockData } from "./index";

/**
 * What a project template carries out of a configured block — the mirror image of
 * `initBlockData`, and the reason the two must be read together: a field added to the contract
 * but not to this function is silently dropped from every template exported afterwards.
 *
 * The fields `BlockParams` leaves out are left out here for the reasons the contract records:
 * the derived label, the settings drawer's own open flag, and view state.
 */
export function deriveTemplateParams(data: BlockData): BlockParams {
  return {
    abundanceRef: data.abundanceRef,
    concentrationColumnRef: data.concentrationColumnRef,
    binColumnRef: data.binColumnRef,
    antigenColumnRef: data.antigenColumnRef,
    sortFractionColumnRef: data.sortFractionColumnRef,
    targetAntigen: data.targetAntigen,
    minReadsPerConcentration: data.minReadsPerConcentration,
    minConcentrationPoints: data.minConcentrationPoints,
    r2ThresholdGood: data.r2ThresholdGood,
    r2ThresholdFailed: data.r2ThresholdFailed,
    nMin: data.nMin,
    nMax: data.nMax,
    hookEffectThresholdBin: data.hookEffectThresholdBin,
    hookEffectThresholdNoBin: data.hookEffectThresholdNoBin,
    hookEffectMinReads: data.hookEffectMinReads,
    customBlockLabel: data.customBlockLabel,
  };
}
