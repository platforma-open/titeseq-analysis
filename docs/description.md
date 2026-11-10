# Overview

The block estimates the apparent binding affinity (Kd,app) of each antibody variant from Tite-Seq data. For each variant, it fits a Hill curve to the signal across the antigen concentrations.

The block has two fit equations:

- **Log fluorescence fit.** The block fits the mean log10 fluorescence of each variant. The Hill exponent is fixed at 1. New blocks use this fit.
- **Mean bin fit.** The block fits the read-weighted mean FACS bin of each variant. The Hill exponent is free. Blocks created before the equation choice existed keep this fit. Without FACS bins, this fit uses the frequency of each variant as the signal.

With the log fluorescence fit, the block gives each variant one affinity class:

- **Good:** Kd,app is inside the titrated range. Its standard error (log10) is 0.5 or less.
- **Partial:** Kd,app is inside the titrated range. Its standard error (log10) is more than 0.5.
- **Below range:** Kd,app is at or below the lowest concentration. The variant binds more tightly than the titration can measure.
- **Above range:** Kd,app is at or above the highest concentration.
- **Failed:** the block cannot fit the variant. The fit failure reason gives the cause.

The mean bin fit uses only Good, Partial and Failed. It rates a single fit by R².

If you select a sample replicate column, the block fits each replicate separately. It then reports the mean of log10 Kd,app across the replicates and the standard error of that mean.

On a nucleotide variant dataset that has a nucleotide to amino-acid linker, for example from Amplicon Profiling, the block fits each nucleotide variant separately. It then reports one Kd,app for each amino-acid variant: the mean of log10 Kd,app across its synonymous nucleotide variants, with the standard error of that mean. The NT Variants page shows the fit of each nucleotide variant.

A mean across replicates or synonymous variants gets its class from the mean, with both fits. It is Good only if its standard error (log10) is 0.5 or less and more than half of the values in the mean were measured inside the titrated range. A value at an edge of the range enters the mean at that edge.

The FACS bin column can hold numbers or labels. With labels, set the bin order in the settings, from the lowest fluorescence to the highest.

The table shows Kd,app in molar units and as -log10 Kd,app. A larger -log10 Kd,app means tighter binding.

The block reads variant abundance from a clonotyping block, for example MiXCR Clonotyping or Amplicon Profiling. Antibody/TCR Lead Selection can use Kd,app and the affinity class as In Vitro preset inputs. You can add a 0 M sample with no antigen. The block uses it to set the Hill baseline, which makes fits of weak binders more reliable.

Kd,app is an apparent affinity, not a thermodynamic Kd. Use it to rank variants in one experiment only. Do not compare it across experiments, display formats (scFv or Fab), or biophysical assays (SPR or BLI). With the mean bin fit, integer bin labels compress the mean-bin scale. This follows the convention in the Tite-Seq literature (Starr et al. 2020; Adams et al. 2016).
