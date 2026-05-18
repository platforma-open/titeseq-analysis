---
'@platforma-open/platforma-open.titeseq-analysis.workflow': patch
'@platforma-open/platforma-open.titeseq-analysis': patch
---

Fix `CIDConflictError` when two TiteSeq blocks with identical settings run in the same project. The Python `fit-curves` binary now sorts every output frame deterministically before writing and pins `maintain_order=True` on every `polars.group_by` whose result reaches a disk write. A prior commit fixed `mean_bin.tsv` only; the same hazard lived in `per_clonotype.tsv`, `fitted_mean_bin.tsv`, and `concentration_value.tsv`. Adds an integration test that runs the binary twice on the same fixture and byte-compares all four outputs.
