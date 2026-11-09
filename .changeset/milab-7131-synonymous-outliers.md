---
"@platforma-open/platforma-open.titeseq-analysis": minor
"@platforma-open/platforma-open.titeseq-analysis.software": minor
"@platforma-open/platforma-open.titeseq-analysis.workflow": minor
"@platforma-open/platforma-open.titeseq-analysis.model": minor
"@platforma-open/platforma-open.titeseq-analysis.ui": minor
---

Flag synonymous nucleotide variants that disagree with the others of their amino-acid variant, and
open a per-variant titration chart from any row.

- Each synonymous variant gets its gap from the mean of the other synonymous variants measured
  inside the range (`pl7.app/vdj/kdLog10SynonymDelta`, log10 Kd), the error of that gap
  (`kdLog10SynonymDeltaSe`), z = gap / error (`kdLog10SynonymZ`) and a flag at |z| >= 2.5
  (`synonymOutlier`), on the nucleotide axis. The amino-acid column `pl7.app/vdj/synonymOutliers`
  names each flagged variant with its direction, gap and z; two disagreeing variants are named as a
  pair, since neither is the odd one.
- The error is sqrt(se^2 + se_others^2), with se_others from the synonymous combination rule. Where
  an SE is missing (`hill_gate_index`, a variant at a wall), it is the dataset sigma,
  1.4826 * MAD of every measured variant's deviation from its amino-acid mean pooled over the run,
  each scaled by sqrt(m / (m - 1)) for the m measured variants in that mean, times
  sqrt(1 + 1 / m_others).
- A variant at a wall enters at the wall value and is tested as a lower bound: its note reads
  "at least X decades, below range" or "above range", without z. It never enters another variant's
  comparison or the dataset sigma.
- The per-nucleotide-variant table is written in every synonymous run, not only when nucleotide
  mutation counts are given.
- The variant panel opens from any row, with or without replicates, and shows the variant's
  titration curves: measured signal per replicate and the fitted curve, pinned to the opened
  variant (all its synonymous variants on the amino-acid grain). The pin is never saved in the
  block's data.
- The NT Variants page switches between all nucleotide variants and synonymous variants only with
  tabs, as the QC page does.
