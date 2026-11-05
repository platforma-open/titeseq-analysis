/** The run-level numbers the software writes for the summary pages (run_stats.json). */
export type RunStats = {
  /** "aa" in synonymous mode, where rows are amino-acid variants; "variant" otherwise. */
  grain: "aa" | "variant";
  variants: number;
  classes: Record<string, number>;
  /** Variants bound at the lowest ("below") and highest ("above") concentration. */
  walls: { below: number; above: number };
  concentrations: { count: number; log10Lowest: number | null; log10Highest: number | null };
  /** Replicate labels, sorted; empty without a replicate column. */
  replicates: string[];
  medianReadsByClass: Record<string, number>;
  /** Parent mode only. */
  parents?: ParentHeadroom[];
};

export type ParentHeadroom = {
  parentId: string;
  log10Kd: number | null;
  log10KdSe: number | null;
  /** Decades from the parent down to the lowest concentration. */
  tightRoom: number | null;
  /** Decades from the parent up to the highest concentration. */
  weakRoom: number | null;
};
