---
'@platforma-open/platforma-open.titeseq-analysis.model': patch
---

Fix `axes sets are disjoint` error in the summary Table when upstream uses `redefine-clonotypes`. The Table builder pulled annotation columns from the result pool with a function predicate that checked axis names only. With redefined clonotypes, two cohorts share the `pl7.app/vdj/clonotypeKey` axis name but have different domains (one carries `pl7.app/redefined-by` + an aminoacid `clonotypeKey/structure`; the other is the original chain-aggregate keyed by nucleotide structure + V/J gene hits). Both passed the predicate and were merged, producing disjoint key spaces. Switched to the declarative anchored-selector form, which matches axes by full spec (name + domain).
