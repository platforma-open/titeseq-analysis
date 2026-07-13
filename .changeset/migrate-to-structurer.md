---
'@platforma-open/platforma-open.titeseq-analysis.model': patch
'@platforma-open/platforma-open.titeseq-analysis.ui': patch
'@platforma-open/platforma-open.titeseq-analysis.workflow': patch
'@platforma-open/platforma-open.titeseq-analysis.software': patch
'@platforma-open/platforma-open.titeseq-analysis': patch
---

Migrate the block onto the `block-tools structure` layout and upgrade the SDK
to latest (model/ui-vue 1.79.x, workflow-tengo 6.x, tengo-builder 4.0.x,
block-tools 2.12.4). Tool-managed config now owns tsconfig, oxlint/oxfmt, turbo,
the block index facade, and per-package deps. The model's block export is
renamed `model` → `platforma` to match the canonical facade contract, and the
test's self-spec import switches to the `from-pack-v2` `TiteseqAnalysisBlockPointer`.
No functional change.