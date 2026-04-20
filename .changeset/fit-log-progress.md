---
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

Stream timestamped progress messages from the Python pipeline to stdout so the Fit Log UI shows live stage transitions (load, validate, normalize, baseline, hook-detect, per-clonotype fit progress, output write). Long fit runs no longer appear stalled — users see elapsed seconds and a running `fitted X/N clonotypes` counter every ~5%.
