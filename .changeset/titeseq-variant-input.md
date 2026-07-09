---
'@platforma-open/platforma-open.titeseq-analysis.model': minor
'@platforma-open/platforma-open.titeseq-analysis.ui': minor
'@platforma-open/platforma-open.titeseq-analysis': minor
---

Accept synthetic-repertoire-profiler datasets as Tite-Seq input

The read-count picker now also offers abundance anchors keyed on
[sampleId, variantKey] (synthetic-repertoire-profiler's pl7.app/readCount),
alongside the existing MiXCR [sampleId, clonotypeKey] / [sampleId, scClonotypeKey]
anchors. The titration-curves plot picks its per-item key axis by excluding the
concentration axis, so facet/grouping defaults resolve for variant datasets too.
The workflow is unchanged — it was already axis-generic.