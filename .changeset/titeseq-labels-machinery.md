---
'@platforma-open/platforma-open.titeseq-analysis.workflow': patch
'@platforma-open/platforma-open.titeseq-analysis.ui': patch
'@platforma-open/platforma-open.titeseq-analysis': patch
---

Fix labels machinery and CID conflicts for multi-instance TiteSeq blocks.

- **Labels:** the workflow trace label now resolves to `customBlockLabel || defaultBlockLabel` (matching clonotype-clustering), so the user's subtitle override propagates to downstream pickers. UI auto-derivation leads with `targetAntigen` (the per-instance differentiator in multi-antigen studies) and drops the bin column — which was identical across instances and pushed the differentiator out of the join. Downstream blocks like Lead Selection now show distinct labels for each instance's "Affinity class" column instead of three identical entries.
- **Determinism:** two TiteSeq blocks with byte-identical settings in the same project no longer hit `CIDConflictError`. The Python `fit-curves` binary now sorts every output frame before writing (`per_clonotype.tsv`, `fitted_mean_bin.tsv`, and a tie-breaker on `concentration_value.tsv`) and pins `maintain_order=True` on every `polars.group_by` whose result reaches a disk write. A prior commit handled `mean_bin.tsv` only; the same hazard lived in the three sibling outputs.
