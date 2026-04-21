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

- Block renamed from "Titeseq Analysis" to "Tite-Seq Analysis" to match the canonical capitalization of the method.
- Friendly input validation: reject all-empty concentration columns early; hide metadata columns that are empty for every sample from the Antigen concentration and FACS bin dropdowns; the bin validator raises a clear "pick a populated column" message when a saved project still references an empty column.
- Dropdown typing: the Antigen concentration picker is filtered to Float/Double metadata columns only — integer-typed bin/replicate columns no longer leak in.
- Dynamic "Target" picker label: the value picker reflects the selected antigen-column's human name (e.g. "Target Sample" when the Antigen Label column is labelled "Sample"). Falls back to "Target antigen" when no column is chosen.
- Informative subtitle: the block subtitle derives from the first three populated inputs joined by " - " (matching Amplicon Alignment), falling back to "Tite-Seq Analysis" when nothing is selected yet; still editable via the Custom label field in Inputs.
- Fit-log progress: the Python pipeline streams timestamped stage transitions to stdout, so the Fit Log UI shows live progress (load, validate, normalize, baseline, hook-detect, per-clonotype fit counter every ~5%, output write) instead of appearing stalled.
- Refreshed block and organization logo assets.
- R14 compliance: `meanBin` and `fittedMeanBin` PColumns carry both the canonical String axis `concentrationStr` (join key, preserving upstream metadata byte-for-byte) and the numeric Long sibling `concentrationAM` (attomolar, drives Graph Maker log-scale rendering). A parametrized unit test guards the `concentrationAM == round(float(concentrationStr) × 1e18)` invariant so the two axes cannot silently desynchronize.
- Block-level integration tests: the `test/` package now exercises Samples-and-Data → import-vdj-data → Tite-Seq Analysis end-to-end against synthetic MiXCR-format fixtures (28 bin-mode samples = 5 clonotypes × 4 FACS bins × 7 concentrations). One consolidated `blockTest` covers option-list population, the six summary PColumn specs, affinityClass distribution (≥1 Good, ≥1 Failed), K_D bounds, and `titrationCurvesPf` composition. CI now runs `pnpm test` for the block; the `workflow/` `test` script carries `--passWithNoTests` since the stub was retired.
