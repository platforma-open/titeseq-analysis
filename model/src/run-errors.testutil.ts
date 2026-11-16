import type { BlockData } from "./index";
import { deriveArgsIssues } from "./index";

/** The reasons Run is disabled, one per line; empty when the block may run. */
export const runErrors = (data: BlockData): string =>
  deriveArgsIssues(data)
    .map((i) => i.message)
    .join("\n");
