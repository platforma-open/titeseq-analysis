---
'@platforma-open/platforma-open.titeseq-analysis.software': minor
'@platforma-open/platforma-open.titeseq-analysis.workflow': minor
'@platforma-open/platforma-open.titeseq-analysis.model': minor
'@platforma-open/platforma-open.titeseq-analysis.ui': minor
'@platforma-open/platforma-open.titeseq-analysis': minor
---

Add optional FACS sort-fraction correction to Mean Bin.

Users can now supply a per-sample sort_fraction metadata column (C_bc/C_c from Adams, Mora, Walczak, Kinney 2016). When present, Mean Bin is sort-yield-corrected and the output carries the annotation pl7.app/titeseq/facsCorrected="true". Absent the column, behaviour is bit-exact with the prior release.

Concentration and Sort-fraction dropdowns disambiguate by data (cross-exclusion + [0,1]-range guard) rather than column names, so sort_fraction won't appear in Antigen Concentration and vice versa.
