# Overview

Estimates apparent binding affinity (Kd,app) per antibody variant from Tite-Seq data. For each clonotype at each antigen concentration, the block computes a signal (mean-bin from per-bin FACS read counts, or clonotype frequency in no-bin mode), fits a Hill equation to the titration curve, and classifies each clonotype as Good, Partial, or Failed.

The block runs downstream of MiXCR Clonotyping, which supplies per-sample, per-clonotype abundance, and upstream of Antibody/TCR Lead Selection, which consumes Kd,app and affinity class as In Vitro preset inputs. An optional 0 M no-antigen control fixes the Hill baseline and improves fit reliability for weak binders.

Kd,app is apparent affinity, not thermodynamic Kd. Values support within-experiment ranking only and are not comparable across experiments, display formats (scFv vs Fab), or biophysical assays (SPR, BLI). Integer bin labels introduce a systematic compression of the mean-bin scale, following the convention in the Tite-Seq literature (Starr et al. 2020; Adams et al. 2016).

# Concentration Units

Numeric concentration appears in two parallel forms. In tables, *Concentration (canonical)* preserves the raw input string (e.g. `5.00E-06`) byte-for-byte. In titration curve plots, the X-axis uses the same data stored as a `Long` in **attomolar** (aM = molar × 10¹⁸); a `Long` axis is required for log-scale rendering without float-roundtrip drift across the workflow's Tengo / Python / TypeScript layers. The two columns reference the same physical concentration.

Conversion: 1 fM = 10³ aM, 1 pM = 10⁶ aM, 1 nM = 10⁹ aM, 1 µM = 10¹² aM, 1 mM = 10¹⁵ aM, 1 M = 10¹⁸ aM. So a curve spanning 10⁶–10¹² aM corresponds to 1 pM – 1 µM.
