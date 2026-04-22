# Overview

Estimates apparent binding affinity (Kd,app) per antibody variant from Tite-Seq data. For each clonotype at each antigen concentration, the block computes a mean-bin signal from per-bin FACS read counts, fits a Hill equation to the titration curve, and classifies each clonotype as Good, Partial, or Failed.

The block runs downstream of MiXCR Clonotyping, which supplies per-sample, per-clonotype abundance, and upstream of Antibody/TCR Lead Selection, which consumes Kd,app and affinity class as In Vitro preset inputs. An optional 0 M no-antigen control fixes the Hill baseline and improves fit reliability for weak binders.

Kd,app is apparent affinity, not thermodynamic Kd. Values support within-experiment ranking only and are not comparable across experiments, display formats (scFv vs Fab), or biophysical assays (SPR, BLI). Integer bin labels introduce a systematic compression of the mean-bin scale, following the convention in the Tite-Seq literature (Starr et al. 2020; Adams et al. 2016).
