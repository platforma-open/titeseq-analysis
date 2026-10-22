---
'@platforma-open/platforma-open.titeseq-analysis.model': patch
'@platforma-open/platforma-open.titeseq-analysis.ui': patch
'@platforma-open/platforma-open.titeseq-analysis.workflow': patch
'@platforma-open/platforma-open.titeseq-analysis.software': patch
'@platforma-open/platforma-open.titeseq-analysis': patch
---

Fix block alerts and settings found in the 2026-10-20 bug hunt

- A failed run shows its cause in the block's alerts, instead of an empty alert list.
- When every clonotype is Failed, Titration Curves says so instead of showing an empty chart.
- The no-bin-mode warning follows the results on screen, not unsaved settings.
- Target antigen is no longer cleared when the block's data is updated from outside the settings panel.
