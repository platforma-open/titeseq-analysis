---
'@platforma-open/platforma-open.titeseq-analysis.workflow': patch
'@platforma-open/platforma-open.titeseq-analysis.ui': patch
'@platforma-open/platforma-open.titeseq-analysis': patch
---

Fix labels machinery for multi-instance TiteSeq blocks. The workflow trace label now resolves to `customBlockLabel || defaultBlockLabel` (matching clonotype-clustering), so the user's subtitle override propagates to downstream pickers. UI auto-derivation leads with `targetAntigen` (the per-instance differentiator in multi-antigen studies) and drops the bin column — which was identical across instances and pushed the differentiator out of the join. Downstream blocks like Lead Selection now show distinct labels for each instance's "Affinity class" column instead of three identical entries.
