---
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

Sort `mean_bin.tsv` rows by (clonotypeKey, concentrationStr) so the output is byte-stable across runs. Polars `group_by` inside `normalize()` returns rows in non-deterministic order, which propagated through `xsv.importFile` to varying Parquet content hashes and triggered CIDConflictError on re-derivation of identical inputs.
