---
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

Address PR #1 review feedback:

- **Parallelize Hill fits** — `_execute_fits` now uses `ProcessPoolExecutor` when survivor count meets a threshold that amortizes pool start-up (default 50). Falls back to the serial loop for small batches. Workflow already allocates 4 CPUs to this step, so large runs with thousands of clonotypes see near-linear speed-ups.
- **Reject non-finite concentrations (R2)** — `validate_concentration_column` now also raises `InputValidationError` on NaN/Inf. Previously these slipped past the `< 0` check and caused opaque downstream failures inside `curve_fit`.
- **Validation warnings go to stdout** — `_validate_inputs` routes WARN lines through the `log()` helper (stdout) instead of stderr so the Tengo workflow's `saveStdoutStream()` surfaces them in the Fit Log UI.
- **Trap `OptimizeWarning` as failure** — `fit_one_clonotype` wraps `curve_fit` in `warnings.catch_warnings()` + `filterwarnings("error", OptimizeWarning)`, so scipy's "covariance could not be estimated" and similar signals produce a `_failure()` result instead of returning a point estimate whose quality scipy itself has flagged as unreliable.
- **Cheaper metadata-uniqueness check (R5)** — `validate_sample_metadata_uniqueness` deduplicates `(sampleId, key_cols)` before grouping. For inputs with millions of rows (clonotypes × samples), the group-by now runs over ≤ n_samples rows instead of the full table.
- **Use `dtype.is_integer()` for bin type check** — replaces the hand-listed `(pl.Int8, ... pl.UInt64)` tuple. `pl.INTEGER_DTYPES` was suggested in review but is deprecated in polars 1.x; `is_integer()` is the current idiomatic API.
