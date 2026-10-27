# @platforma-open/platforma-open.titeseq-analysis.kind

## 1.1.0

### Minor Changes

- 03ca874: Update the SDK and accept a subset filter (e.g. a Repertoire Labeling label or Lead Selection picks) with the input read counts. The dataset list is built with buildDatasetOptions (requireEnrichments: false). The subset restricts the input before the fit, so the per-bin depths and clonotype frequencies are computed within the subset, and Kd,app values can differ from a full-data run. Outputs of a subset run carry pl7.app/inputSubset in their domain and a trace rooted in the subset, so their labels name it; the default subtitle names the subset too. Unfiltered runs are unchanged.
