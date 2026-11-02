import type { Equation } from "@platforma-open/platforma-open.titeseq-analysis.model";
import { useApp } from "../app";

/** Writing fit settings into block data. */
export function useFitSettings() {
  const app = useApp();

  /**
   * Writes the equation and the Affinity vs Fit Quality y-axis scale. The y column is the SE under
   * log_fluorescence (log scale) and hillPlotPosition under hill_gate_index (linear: it holds -1).
   */
  function setEquation(equation: Equation): void {
    app.model.data.equation = equation;
    const graph = app.model.data.graphStateAffinityVsFit;
    graph.axesSettings = {
      ...graph.axesSettings,
      axisY: {
        ...graph.axesSettings?.axisY,
        scale: equation === "log_fluorescence" ? "log" : "linear",
      },
    };
  }

  return { setEquation };
}
