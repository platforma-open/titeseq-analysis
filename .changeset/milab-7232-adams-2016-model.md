---
'@platforma-open/platforma-open.titeseq-analysis.software': minor
---

Fit bin-mode K_D with the Adams et al. 2016 (eLife 5:e23156) estimator by default

The new default model, `poisson_log`, works on read counts:
- Poisson likelihood on reads.
- Mean log fluorescence modelled as log10(B + A·c/(c+K)), with Hill n = 1.
- Bin probabilities weighted by the cells sorted into each bin, and reads rescaled where they exceed cells.
- One global 0 M baseline.
- A minimum-plateau amplitude prior, so flat curves read as K_D above the range.

It reproduces the paper:
- Per-variant K_D matches the authors' published values at r = 0.95–0.99.
- WT averages 10^-8.85 (paper: 10^-8.87).
- Table I vs flow cytometry: r² = 0.88 (paper: 0.89).

New inputs:
- `--cell-count-column` (cells sorted per sample).
- `bin_width_decades`, `min_plateau_fraction`, `baseline_fluorescence`, `allow_uncorrected` and `kd_sd_good`.

Bin mode now requires cell counts or sort fractions unless `allow_uncorrected` is set. It also requires a 0 M control unless `baseline_fluorescence` is given.

New outputs:
- `kdLog10Sd`: posterior SD of log10 K_D.
- `kdCensor`: `below` or `above` when K_D is a bound at the titration limit.

K_D is reported within the titration range. The previous linear Hill model with free n stays available as `fit_model: hill_free_n`.

The linear Hill model (hill_free_n) now reports a flat signal as Failed (convergence_failure) instead of Partial with an out-of-range K_D.
