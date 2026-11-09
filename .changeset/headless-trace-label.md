---
"@platforma-open/platforma-open.titeseq-analysis.model": patch
"@platforma-open/platforma-open.titeseq-analysis.ui": patch
"@platforma-open/platforma-open.titeseq-analysis.workflow": patch
---

Compute the block label without the UI

The trace label is now the custom block label or "Tite-Seq Analysis (<target antigen>)",
computed by the model, so opening the block no longer changes its arguments and a block
that was never opened gets the same label. The abundance column's name moves to the
subtitle only, written when the user picks the abundance; it was already part of the
inherited upstream trace. Existing blocks show as stale once after the update because
the arguments change shape.
