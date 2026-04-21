---
'@platforma-open/platforma-open.titeseq-analysis': major
'@platforma-open/platforma-open.titeseq-analysis.model': major
'@platforma-open/platforma-open.titeseq-analysis.software': major
'@platforma-open/platforma-open.titeseq-analysis.test': major
'@platforma-open/platforma-open.titeseq-analysis.ui': major
'@platforma-open/platforma-open.titeseq-analysis.workflow': major
---

Initial public release of the Tite-Seq Analysis block.

Fits Hill curves per clonotype against MiXCR abundance across concentrations, emits Kd,app with a confidence class (Good / Failed), and renders Titration Curves, Kd Distribution, and Affinity vs Fit Quality tabs.

- Spec compliance: R5 (sample-metadata uniqueness), R9b (top-1/top-2 hook gate), R10 (Hill baseline bounds), R14 (dual concentration axes), R15 (Titration Curves default layout with affinity filter), R17 (Failed clonotypes park at Kd = -1.0 and render on Affinity vs Fit Quality).
- R14 invariant: `meanBin`/`fittedMeanBin` carry both `concentrationStr` (String join key) and `concentrationAM` (Long, attomolar, log-scale axis); a parametrized test guards the round-trip.
- Inputs: reject empty/NaN/Inf concentration columns; hide all-empty metadata columns from pickers; restrict Antigen concentration to Float/Double; Target picker label follows the selected Antigen column.
- UX: block subtitle derives from the first three populated inputs; each tab shows a single page title; Langmuir tooltip on the Hill coefficient PColumn; input labels and warnings repositioned for clarity.
- Fit Log streams timestamped stage progress (load → validate → normalize → baseline → hook-detect → per-clonotype fit → write) and surfaces validator warnings.
- Performance: Hill fits parallelize above 50 clonotypes via `ProcessPoolExecutor`; window functions replace self-joins in normalization; single-pass boundary counters in validators; metadata-uniqueness check scales with n_samples.
- Robustness: scipy `OptimizeWarning` marks fits Failed; attomolar ceiling guarded; concentration cast defensively to Float64.
- Architecture: migrated to `BlockModelV3`; concentration/bin/antigen columns anchored; Python runtime pinned to LTS-CPU wheels.
- Tests: block-level integration exercises Samples-and-Data → import-vdj → Tite-Seq against a 28-sample synthetic fixture (5 clonotypes × 4 bins × 7 concentrations); a deterministic e2e corpus suite bounds numeric drift on per-clonotype columns.
