---
"@platforma-open/platforma-open.titeseq-analysis": patch
"@platforma-open/platforma-open.titeseq-analysis.model": patch
"@platforma-open/platforma-open.titeseq-analysis.ui": patch
---

- Update graph-maker from 1.3.1 to 1.10.5. A saved chart now takes a default filter that lost an
  entry: a project last run under the mean bin fit no longer keeps the Failed-only filter on the
  Affinity vs Fit Quality page, which left that plot empty under the log fluorescence fit. The
  Kd Distribution tooltip no longer shows `Symbol(noGrouped)`.
- Settings shows on the main page only. The other pages keep Logs.
