---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.kind": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
---

Update the SDK and accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) with the input read counts. The dataset list is built with buildDatasetOptions (requireEnrichments: false). The subset restricts the input before the fit, so the per-bin depths and clonotype frequencies are computed within the subset, and Kd,app values can differ from a full-data run. Outputs of a subset run carry pl7.app/inputSubset in their domain and a trace rooted in the subset, so their labels name it; the default subtitle names the subset too. Unfiltered results are unchanged, but existing blocks re-run once after the update: the workflow template gets a new hash override.
