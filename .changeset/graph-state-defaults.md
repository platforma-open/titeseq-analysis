---
'@platforma-open/platforma-open.titeseq-analysis.model': minor
'@platforma-open/platforma-open.titeseq-analysis': minor
---

Persist in-use plot configuration as block defaults for Kd Distribution and Affinity vs Fit Quality.

Kd Distribution opens with a log y-axis and green bin fill; Affinity vs Fit Quality opens with a log x-axis, filtered to Failed clonotypes with Hill coefficient ≥ 0, and coloured by fit-failure reason (low_r2 → green, n_out_of_range → purple). Values bind to the stable summaryPf output columns so the defaults transfer to any new block or project. Titration Curves was already aligned with the intended default.
