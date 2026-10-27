---
'@platforma-open/platforma-open.titeseq-analysis.software': minor
'@platforma-open/platforma-open.titeseq-analysis.workflow': minor
'@platforma-open/platforma-open.titeseq-analysis.model': minor
'@platforma-open/platforma-open.titeseq-analysis.kind': minor
'@platforma-open/platforma-open.titeseq-analysis.ui': minor
---

Report every converged fit, route the curve model on the inputs, and never refuse to run

The block no longer discards a fit it managed to make, and no longer stops because an input is
missing. Each extra column the user supplies now buys a better result instead of unlocking a run.

**Classification.** A converged fit is always reported with its class. The low-R² and
out-of-range-Hill-coefficient Failed verdicts are gone; those fits are Partial. Restoring that gate
was measured to discard 215, 209 and 217 of the 250 tightest binders on the three Adams 2016
replicates. The hook-effect detector is off by default and the Hill exponent is fixed at 1, which
is what makes the detector redundant.

**Routing.** The curve model follows the bound columns: sorted cell counts (or sort fractions)
select the Adams 2016 measurement model, their absence the ranking fit. Feeding uncorrected bins to
the measurement model is measurably worse than ranking them, so more data has to switch the model
rather than only refine it. `fitModelOverride` forces a model by hand. A project saved before this
release is routed like a new one.

**Bounds.** A K_D outside the titrated concentrations is clamped to the limit and flagged, never
failed. Below the range it carries its own class, `CensoredBelow`, outside the default
Good + Partial cutoff, because the clamped value is the smallest in the table and would otherwise
head a downstream ranking as an unorderable block of bounds.

**The 0 M baseline** is computed whenever a 0 M sample exists, from the clonotypes deep enough to
carry it rather than from every clonotype equally. The old estimator converged on a wrong number as
shallow sequencing-error clonotypes came to dominate. The new one lands within 1.7%, 2.7% and 2.4%
of the authors' hand values and warns when the 0 M set looks contaminated.

**Two optional FACS inputs.** `binFluorescenceColumnRef` supplies the measured mean log10
fluorescence of each sample's gate, replacing one fixed bin scale with one value per bin per
concentration; a gate's real mean moves with concentration, and a fixed value reads every midpoint
tighter than it is. `gateWidthColumnRef` supplies the bin width. Either one, if unusable, warns and
falls back rather than failing.

**Removed:** the two hard stops (Adams model without cell counts, computed baseline without opt-in),
the `allowUncorrected` setting and the `--allow-computed-baseline` flag. In their place the block
warns about what it did and what one more input would buy.

Measured on the three Adams 2016 replicates with bins and a 0 M control but no cell counts, scored
on the rows the user is shown: 330 / 0.888 / 15 becomes 1,706 / 0.909 / 236 on replicate 1,
97 / 0.822 / 1 becomes 1,698 / 0.889 / 218 on replicate 2, and 467 / 0.878 / 23 becomes
1,388 / 0.946 / 243 on replicate 3 (rows shown / rank agreement against the authors' values /
retained of their 250 tightest).
