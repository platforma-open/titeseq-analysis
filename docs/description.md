# Overview

The Titeseq Analysis block reconstructs apparent binding affinity (K_D,app) for antibody variants from titeseq sequencing data. For each clonotype at each antigen concentration, the block computes a mean-bin signal from per-bin read counts, fits a Hill equation to the resulting titration curve, and ranks variants by measured binding affinity.

The block sits in the antibody-discovery pipeline downstream of MiXCR clonotyping — which supplies per-sample, per-clonotype abundance — and upstream of Lead Selection, which consumes K_D,app and affinity class as In Vitro preset inputs.

# Inputs

- **Abundance** — MiXCR-derived per-sample, per-clonotype read counts (PColumn with axes `[sampleId][clonotypeKey]`).
- **Concentration Metadata** — per-sample antigen concentration; the column label (e.g. `nM`, `µM`) is propagated to the K_D,app unit annotation.
- **Bin Metadata** (optional) — per-sample integer FACS bin index. When provided, the block runs in bin mode and the signal is mean bin; when absent, the block runs in no-bin mode and the signal is clonotype frequency (not comparable to bin-derived K_D,app).
- **Antigen Metadata** (optional) — per-sample antigen identifier. When provided, the block filters to a single user-selected `targetAntigen` before fitting, so one titeseq run covering multiple antigens can be analyzed one antigen at a time.
- **0 M No-Antigen Control** (recommended) — fixes a global baseline and reduces the Hill fit from four to three free parameters, improving reliability for weak binders.

# Outputs

Per clonotype:

- **K_D,app** — apparent dissociation constant.
- **Hill Coefficient n** and **Curve Fit R²** — fit diagnostics.
- **Affinity Class** — Good / Partial / Failed. Primary hard filter for Lead Selection.
- **Fit Failure Reason** — machine-readable reason for Failed clonotypes (`insufficient_reads`, `insufficient_points`, `non_monotonic_signal`, `convergence_failure`, `low_r2`, `n_out_of_range`).
- **K_D Out-of-Range Flag** — true when K_D,app falls outside the experimental concentration grid (extrapolation warning, not suppression).

Per clonotype and concentration:

- **Mean Bin** — observed signal, or clonotype frequency in no-bin mode.
- **Fitted Mean Bin** — Hill curve evaluated at each experimental concentration point, used as the overlay in the titration-curve plot.

# Visualizations

- **Titration Curves** — mean-bin dots with the fitted Hill sigmoid overlaid per clonotype; concentration axis on log scale.
- **K_D Distribution** — histogram of log₁₀ K_D,app across Good and Partial clonotypes.
- **Affinity vs Fit Quality** — scatter of K_D,app (log scale) against Hill coefficient n, colored by affinity class. Failed clonotypes are included so users can tune `nMin`, `nMax`, and R² thresholds without rerunning the workflow.
- **Table** — all clonotypes with all output columns, sortable.

# Caveats

K_D,app is apparent affinity, not thermodynamic K_D. Values are valid for within-experiment ranking and should not be compared across experiments, display formats (scFv vs Fab), or against SPR or BLI biophysical assays. Integer bin labels introduce a systematic compression of the mean-bin scale that is accepted in the titeseq literature (Starr et al. 2020; Adams et al. 2016).
