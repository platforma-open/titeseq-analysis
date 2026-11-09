import type { BlockData } from "@platforma-open/platforma-open.titeseq-analysis.model";
import {
  EXPLORER_CURVES_DEFAULT_STATE,
  EXPLORER_SCATTER_DEFAULT_STATE,
} from "@platforma-open/platforma-open.titeseq-analysis.model";

/**
 * Starts both Explorer curve charts again from their defaults. They share the class pick, which
 * reaches each as its default filter, and a chart applies its defaults once.
 */
export function resetExplorerCurves(
  data: Pick<BlockData, "graphStateExplorerCurves" | "graphStateExplorerScatter">,
) {
  data.graphStateExplorerCurves = { ...EXPLORER_CURVES_DEFAULT_STATE };
  data.graphStateExplorerScatter = { ...EXPLORER_SCATTER_DEFAULT_STATE };
}
