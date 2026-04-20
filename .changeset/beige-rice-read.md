---
'@platforma-open/platforma-open.titeseq-analysis': major
'@platforma-open/platforma-open.titeseq-analysis.model': major
'@platforma-open/platforma-open.titeseq-analysis.software': major
'@platforma-open/platforma-open.titeseq-analysis.test': major
'@platforma-open/platforma-open.titeseq-analysis.ui': major
'@platforma-open/platforma-open.titeseq-analysis.workflow': major
---

Initial public release of the Tite-Seq Analysis block.

Highlights:

- R14 compliance: `meanBin` and `fittedMeanBin` PColumns carry both the canonical String axis `concentrationStr` (join key, preserving upstream metadata byte-for-byte) and the numeric Long sibling `concentrationAM` (attomolar, drives Graph Maker log-scale rendering).
- Friendly input validation: reject all-empty concentration columns early; hide metadata columns that are empty for every sample from the Antigen concentration and FACS bin dropdowns; the bin validator raises a clear "pick a populated column" message when a saved project still references an empty column.
- Dropdown typing: the Antigen concentration picker is filtered to Float/Double metadata columns only — integer-typed bin/replicate columns no longer leak in.
- Dynamic "Target" picker label: the value picker reflects the selected antigen-column's human name (e.g. "Target Sample" when the Antigen Label column is labelled "Sample"). Falls back to "Target antigen" when no column is chosen.
- Informative auto-subtitle: the block subtitle derives from the selected read-count column and mode (`<dataset> · bin mode` / `<dataset> · frequency mode`), editable via the Custom label field in Inputs.
- Fit-log progress: the Python pipeline streams timestamped stage transitions to stdout, so the Fit Log UI shows live progress (load, validate, normalize, baseline, hook-detect, per-clonotype fit counter every ~5%, output write) instead of appearing stalled.
- Refreshed block and organization logo assets.
