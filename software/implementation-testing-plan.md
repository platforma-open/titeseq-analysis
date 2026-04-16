# Affinity Profiling — Python implementation & testing plan

**Scope:** the Python fitting script that runs inside the affinity-profiling block. The Tengo workflow hands this script a normalized reads table; the script emits per-clonotype K_D,app / n / R²_w / class / failure reason / kdOutOfRange and per-(clonotype, concentration) meanBin / fittedMeanBin values. Everything outside the script boundary (Tengo workflow, model, UI) is out of scope for this plan.

**Reference prior art:** `blocks/antibody-sequence-liabilities/liabilities-calc-script` — our canonical "polars + uv + pyproject.toml + pytest" pattern. We follow that layout.

**Authoritative specs:**
- `./README.md` (requirements R1–R21, defaults, failure modes, rationale)
- `./pcolumn-spec.md` (output column names, value types, annotations, failure-reason enum)

---

## 1 — Package layout

```
affinity-profiling-calc-script/
├── pyproject.toml                    # uv-managed; python >= 3.12
├── uv.lock
├── package.json                      # Platforma block-software manifest
├── src/
│   ├── main.py                       # argparse entrypoint; orchestrates pipeline
│   ├── requirements.txt              # runtime deps (pip toolset for pl-pkg)
│   ├── io_layer.py                   # read inputs (parquet/tsv) + canonical concentration key handling
│   ├── normalization.py              # R7, R7b: freq + mean bin / frequency signal (polars-native)
│   ├── floor_filter.py               # R8, R9: weight computation + read-floor + insufficient-* classification
│   ├── hook_effect.py                # R9b: pre-fit non-monotonic signal check
│   ├── hill_fit.py                   # R10, R11: reparametrized Hill model + weighted NLS + weighted R²
│   ├── baseline.py                   # R6: global B from c=0 mean-bin (or no-bin signal)
│   ├── classify.py                   # R12: R²_w × n → Good/Partial/Failed + fitFailureReason
│   ├── output_build.py               # R13, R14, R14b: assemble output frames with axis canonicalization
│   └── constants.py                  # default params, DELTA, FAILURE_REASONS, annotation keys
└── tests/
    ├── conftest.py                   # shared fixtures: synthetic reads builders, param dicts
    ├── data/                         # small fixed TSVs for golden-value tests (pinned Hill outputs)
    ├── unit/
    │   ├── test_normalization.py
    │   ├── test_floor_filter.py
    │   ├── test_hook_effect.py
    │   ├── test_baseline.py
    │   ├── test_hill_fit.py
    │   ├── test_classify.py
    │   ├── test_output_build.py
    │   └── test_io_layer.py
    ├── integration/
    │   ├── test_bin_mode_pipeline.py     # end-to-end bin mode via main.run()
    │   ├── test_no_bin_mode_pipeline.py
    │   └── test_antigen_filter.py
    └── regression/
        └── test_synthetic_titeseq.py     # Poisson-simulated reads; spec R-test: ≥90% K_D within 10%
```

**Why module-per-requirement-cluster:** each file owns one of R6/R7/R8/R9/R9b/R10/R11/R12/R13 so a test failure points at one file. Keeps the hot path (`normalization` + `hill_fit`) narrow and independently profilable.

---

## 2 — Tooling (uv + pyproject.toml)

Mirrors `liabilities-calc-script/pyproject.toml` with the fitting stack added.

```toml
[project]
name = "affinity-profiling-calc-script"
version = "0.0.0"
requires-python = ">=3.12"

[dependency-groups]
dev = [
    "polars>=1.39.0",
    "numpy>=2.0.0",
    "scipy>=1.14.0",
    "pyarrow>=18.0.0",
    "pytest>=9.0.2",
    "pytest-cov>=6.0.0",
]

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src"]
addopts = "-ra --strict-markers"

[tool.ruff]
line-length = 120
target-version = "py312"
# (rest matches liabilities-calc-script)
```

Runtime `src/requirements.txt` (pip, consumed by `pl-pkg build`) pins the same versions (polars, numpy, scipy, pyarrow) — dev extras like pytest are excluded.

**Local commands:**
```
uv sync                              # install deps
uv run pytest tests/                 # run full suite
uv run pytest tests/unit -v          # unit only (fast)
uv run pytest --cov=. tests/         # coverage
uv run pytest -k hill_fit            # single module
uv run pytest -m slow                # regression (see markers below)
```

**Markers:**
- `@pytest.mark.slow` — regression suite with Poisson-sim data (~30 s; excluded from default run via `-m "not slow"` in CI fast job).
- `@pytest.mark.hot_path` — perf smoke tests (see §7).

---

## 3 — Implementation principles

**Polars-native everywhere except the fitting kernel.** All per-clonotype reductions (freq, mean bin, w_j sums, c=0 filter, canonical axis join) are `pl.Expr` pipelines — no `map_elements`, no per-row Python, no `.to_pandas()` round-trips. scipy's curve_fit is the only necessary Python loop.

**Scalar inputs/outputs at the module boundary.** The fitting kernel `fit_one_clonotype(signal, concs, weights, baseline_fixed) -> FitResult` is a pure function over numpy arrays. This lets unit tests pin exact numerical outputs. The per-clonotype loop lives in `hill_fit.py::fit_all()` which extracts arrays from a polars groupby and is the only hot loop.

**Never mutate input frames.** All transforms return new frames. Enables parameter-sweep tests without resetting fixtures.

**Fail fast on invalid inputs.** `io_layer.validate_concentration_column` raises `ValueError` on negative non-zero values (R2) before any fitting starts. Matches spec's "user-facing error" wording.

**Concentration axis canonicalization (R14).** Concentrations flow through as `{canonical_str: str, value: float}` pairs. `output_build` joins fitted values on `canonical_str` so float-serialization drift cannot break the join. This is the single subtlest requirement and has dedicated tests (§6.9).

**No silent clamps.** R²_w can be negative — we let it flow through (spec explicitly forbids clamping to 0).

---

## 4 — Hot paths & performance

Realistic workload: ~1k–10k clonotypes × 5–8 concentrations × 1–8 bins. Targets:

| Path | Target | Technique |
|---|---|---|
| Mean bin / frequency (R7, R7b) | O(rows) single polars pass | `groupby([clonotype, conc]).agg(...)` with expressions; no Python UDFs |
| Weighted R² (R11) | Vectorized polars/numpy; computed post-fit on arrays | `np.average(..., weights=w)` — no Python loop |
| Hill fit (R10) | ~1–10 ms per clonotype × N clonotypes; dominates wall time | scipy `curve_fit` with `method="trf"`, bounded; arrays pre-extracted per group |
| Output assembly (R13, R14) | Single polars concat + `pivot` | No per-clonotype DataFrame allocations |
| Axis canonicalization (R14) | String-keyed left join | `pl.DataFrame` with explicit schema |

**Perf guardrails in the test suite** (see §6.10): one test asserts 1000-clonotype × 8-conc × 4-bin fitting finishes under a generous wall-clock budget (e.g. 30 s on CI). Not a benchmark — just a regression alarm.

**Conscious trade-off:** scipy does not release the GIL inside curve_fit, so we do **not** prematurely parallelize. If profiling shows fitting is the bottleneck, a later optimization can use `multiprocessing.Pool` over clonotypes — but only after measurement. The plan does not pre-commit to parallelism.

---

## 5 — Requirement → module → test mapping

Every requirement from README.md § Requirements has at least one positive and (where meaningful) one negative test.

| Req | Module | Positive test(s) | Negative / edge test(s) |
|---|---|---|---|
| R1 abundance input shape | `io_layer` | `test_reads_wide_long_shape_ok` | `test_missing_clonotype_col_raises`, `test_missing_sample_col_raises` |
| R2 concentration validation | `io_layer` | `test_positive_floats_accepted`, `test_zero_with_bin_accepted` | `test_negative_concentration_raises`, `test_zero_without_bin_assignment_warns` (the ambiguous-0 M-control case — spec R2 explicitly requires a warning) |
| R3 optional bin column | `io_layer` | `test_bin_mode_detected_when_col_present`, `test_no_bin_mode_when_absent`, `test_non_consecutive_bin_labels_ok` (labels `[1,2,5,8]` — spec: "does not assume a fixed number of bins"; `max_bin_label` must be the actual max, not the bin count) | `test_bin_non_integer_raises`, `test_bin_negative_raises`, `test_bin_zero_raises` (spec: positive integer) |
| R4 optional antigen column | `io_layer` | `test_antigen_filter_keeps_target`, `test_antigen_filter_drops_others`, `test_targetAntigen_set_without_antigenColumnRef_warns` (spec: warn, not error) | `test_antigenColumnRef_without_targetAntigen_raises` (spec: model-level user-facing error), `test_targetAntigen_not_in_column_raises` |
| R5 sample ↔ metadata uniqueness | `io_layer` | `test_each_sampleId_one_concentration_ok`, `test_each_sampleId_one_bin_ok`, `test_each_sampleId_one_antigen_ok` | `test_sampleId_with_two_concentrations_raises`, `test_sampleId_with_two_bins_raises`, `test_sampleId_with_two_antigens_raises` |
| R6 global baseline B from c=0 | `baseline` | `test_B_equals_arithmetic_mean_of_c0_meanbins`, `test_B_ignores_filtered_c0_clonotypes` | `test_no_c0_returns_None_B`, `test_all_c0_filtered_returns_None_B`, `test_c0_meanbin_exactly_zero_counted_in_B` (zero is a valid observation, not a sentinel) |
| R7 mean bin formula | `normalization` | Pinned values: §6.1 | `test_zero_denominator_sample_excluded`, `test_single_bin_equals_that_bin_label` |
| R7b no-bin signal | `normalization` | `test_no_bin_freq_equals_reads_over_total` | `test_zero_total_reads_at_conc_excluded` |
| R8 weight + floor filter | `floor_filter` | `test_w_j_equals_bin_read_sum`, `test_floor_excludes_points_below_threshold` | `test_exactly_at_floor_included`, `test_weight_zero_point_dropped` |
| R9 insufficient_reads / insufficient_points | `floor_filter` | `test_zero_points_remaining_marks_insufficient_reads`, `test_under_min_points_marks_insufficient_points` | `test_exactly_min_points_not_marked_insufficient` |
| R9b hook effect | `hook_effect` | Pinned values: §6.3 | `test_below_hook_min_reads_skipped`, `test_monotonic_not_flagged`, `test_small_drop_not_flagged` |
| R10 Hill fit (reparam + bounds) | `hill_fit` | Pinned synthetic fits: §6.5 | `test_flat_curve_yields_convergence_failure`, `test_top_minus_baseline_below_delta_yields_failure`, `test_bounds_enforced` |
| R11 weighted R² | `hill_fit` | Pinned values: §6.6 | `test_negative_R2_not_clamped`, `test_zero_variance_y_returns_nan_sentinel` |
| R12 classification | `classify` | Parametrized truth table: §6.7 | `test_boundary_values_are_inclusive_exclusive_per_spec` |
| R13 output PColumn set | `output_build` | `test_all_required_columns_present`, `test_axes_correct_per_column` | `test_failed_clonotype_has_null_kd`, `test_failed_clonotype_has_failure_reason` |
| R14 meanBin / fittedMeanBin | `output_build` | `test_meanBin_axes_clonotype_concentration`, `test_fittedMeanBin_at_experimental_concs_only` | `test_c0_excluded_from_output`, `test_failed_fit_fittedMeanBin_null` |
| R14b kdOutOfRange flag | `output_build` | `test_kd_in_range_flag_false`, `test_kd_above_max_flag_true`, `test_kd_below_min_flag_true` | `test_kd_exactly_at_min_is_in_range`, `test_kd_exactly_at_max_is_in_range`, `test_failed_fit_kdOutOfRange_null` |
| R15–R19b UI | — | (out of scope for Python; UI tests live in `ui/` Vitest suite) | |
| R20, R21 integration | — | (block-level integration test — workflow/test/src/wf.test.ts) | |

"All clonotypes Failed" dataset-level edge case: `integration/test_bin_mode_pipeline.py::test_all_failed_completes_without_error`.

---

## 6 — Test details: pinned values & parametrization

Formula tests use three families of inputs: **pinned analytical** (hand-computed expected value), **zero/ceiling boundary** (smallest/largest defensible input), **parametrized truth table** (all branches of a decision). Every formula test comment states the invariant it guards.

### 6.1 Mean bin formula (R7) — pinned

```python
# Two clonotypes, one concentration, 4 bins; hand-computed expected mean bin.
# Guards against accidentally switching to raw-count formula (Adams/Starr) — spec deliberately differs.
@pytest.mark.parametrize("reads_per_bin, depth_per_bin, expected_mean_bin", [
    # clonotype A: 10/1000, 20/500, 5/200, 5/100 → freqs [0.010, 0.040, 0.025, 0.050]
    #   mean_bin = (1·0.01 + 2·0.04 + 3·0.025 + 4·0.05) / (0.01+0.04+0.025+0.05)
    #           = 0.365 / 0.125 = 2.92
    ([10, 20, 5, 5], [1000, 500, 200, 100], pytest.approx(2.92, abs=1e-9)),
    # Single-bin clonotype: only bin 3 has reads → mean_bin ≡ 3.0 (pass-through)
    ([0, 0, 7, 0], [1000, 500, 200, 100], 3.0),
    # Uniform frequency across 4 bins → mean_bin = 2.5 (arithmetic centre)
    ([100, 50, 20, 10], [100, 50, 20, 10], 2.5),
])
def test_mean_bin_pinned(reads_per_bin, depth_per_bin, expected_mean_bin): ...
```

Boundary cases:
- **Zero bin reads for one bin** — that bin contributes 0 to numerator and 0 to denominator (frequency = 0); other bins dominate.
- **All-zero reads for a clonotype at a concentration** — excluded by R8 floor, so mean_bin never computed (asserted in `floor_filter` test).
- **Ceiling** — max_bin_label = 8 (spec allows up to 8 bins); a clonotype with reads only in bin 8 yields mean_bin = 8.0.

### 6.2 No-bin frequency signal (R7b) — pinned

```python
@pytest.mark.parametrize("reads_clonotype, reads_total, expected_freq", [
    (100, 10_000, 0.01),
    (1, 10_000, 0.0001),          # floor-threshold-relevant low end
    (10_000, 10_000, 1.0),        # ceiling: every read is this clonotype
    (0, 10_000, 0.0),             # zero numerator
])
def test_no_bin_signal_pinned(reads_clonotype, reads_total, expected_freq): ...
```

Separate test asserts `ZeroDivisionError` path: total_reads_at_conc = 0 → point is R8-floored before signal is computed; signal function never sees it. Test invokes `compute_no_bin_signal(reads_clonotype=5, reads_total=0)` and expects `pl.Null` (or filtered-out row), not an exception.

### 6.3 Hook-effect detection (R9b) — pinned

```python
# hookEffectThresholdBin default = 0.2; hookEffectMinReads default = 20.
# Invariant: signal at top-1 conc must be >= signal at top-2 minus threshold.
@pytest.mark.parametrize("top2_signal, top1_signal, top2_reads, top1_reads, expected_flag", [
    (3.0, 2.5, 100, 100,   False),  # drop 0.5 > 0.2 threshold — flag in bin mode
    (3.0, 2.5, 100, 100,   True ),  # (same data, bin mode) → param'd separately; see mode axis
    (3.0, 2.85, 100, 100,  False),  # drop 0.15 < 0.2 — not flagged
    (3.0, 3.2,  100, 100,  False),  # signal rises — not flagged
    (3.0, 2.0, 100, 10,    False),  # top1_reads < hookEffectMinReads=20 → skipped
    (3.0, 2.0, 10,  100,   False),  # top2_reads < 20 → skipped
    (3.0, 3.0, 100, 100,   False),  # flat (zero drop) — not flagged (strict >)
])
def test_hook_effect_bin_mode(...): ...

# Mirror case for no-bin mode with hookEffectThresholdNoBin=0.02; separate test function for clarity.
```

Boundary: drop exactly equal to threshold (0.2) — spec wording "signal drop at max concentration triggering non_monotonic_signal" is ambiguous on inclusivity. Plan: use strict `>` (drop > threshold flags), documented in a comment referencing this decision. Open question flagged in §8.

### 6.4 Global baseline B (R6) — pinned

```python
# Invariant: B is the arithmetic mean of per-clonotype mean_bin at c=0, ignoring clonotypes whose c=0 was floored.
def test_B_pinned():
    # 3 clonotypes: mean_bin at c=0 = [1.2, 1.4, 1.6] → B = 1.4
    # 4th clonotype's c=0 point below read floor → excluded → B still = 1.4 (not 1.4*3/4)
    ...

# Edge: no c=0 sample at all → B=None → fit uses 4-param (free baseline).
def test_B_none_when_no_c0():
    ...

# Edge: c=0 samples exist but every clonotype's c=0 was floored → B=None.
def test_B_none_when_all_c0_filtered():
    ...
```

### 6.5 Hill fit kernel (R10) — pinned + round-trip

The Hill function with B fixed:
`y(x) = B + exp(amplitude) × xⁿ / (K_D^n + xⁿ)`

**Round-trip test** (strongest guarantee of correctness): generate noiseless y values from `(B=1.0, amplitude=log(3.0), K_D=10.0, n=1.0)` at concentrations `[0.1, 1, 3, 10, 30, 100, 300, 1000]`, fit, assert recovered params within 1e-6. Repeat for `n=2.0` and `n=0.5` (cooperativity range).

```python
@pytest.mark.parametrize("true_kd, true_n, true_amp, conc, expected_kd_abs_err", [
    (10.0, 1.0, math.log(3.0), LOG_CONCS_5p, 1e-6),    # canonical
    (10.0, 2.0, math.log(3.0), LOG_CONCS_5p, 1e-6),    # Hill n=2 cooperativity
    (10.0, 0.5, math.log(3.0), LOG_CONCS_5p, 1e-6),    # n=0.5 anti-cooperativity
    (0.1,  1.0, math.log(3.0), LOG_CONCS_5p, 1e-4),    # low-K_D limit (bound: kd > 0)
    (1000.0, 1.0, math.log(3.0), LOG_CONCS_5p, 1e-4),  # high-K_D limit
])
def test_hill_fit_noiseless_roundtrip(true_kd, true_n, true_amp, conc, expected_kd_abs_err):
    # Asserts recovered K_D within tolerance of the true K_D used to generate the signal.
    # Guards against: wrong reparametrization, sign errors in amplitude, wrong Hill algebra.
    ...
```

**Flat curve (negative case):** generate y = constant + tiny jitter → assert `FitResult.failed == True` and `reason == "convergence_failure"`. This specifically exercises the `top - baseline < δ` path (R10 spec; δ = 0.05 per spec).

**Bounds enforcement:** inject a `kd = -1` initial guess → verify kernel rejects/clips to positive and fit still succeeds (because bounds force kd > 0). Injects `n = 100` initial guess → verify fit converges within bounds `n ∈ [0.1, 10]`.

**Amplitude upper bound depends on mode:** bin mode uses `log(max_bin_label − 1)`; no-bin uses `log(0.95)`. Parametrize over mode.

### 6.6 Weighted R² (R11) — pinned

```python
# Invariant: R²_w = 1 − Σ(w·(y−ŷ)²) / Σ(w·(y−ȳ_w)²), ȳ_w = Σ(w·y)/Σ(w).
# Hand-computed expected values.
@pytest.mark.parametrize("y, y_hat, w, expected_r2_w", [
    ([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 1.0, 1.0], 1.0),     # perfect fit
    ([1.0, 2.0, 3.0], [2.0, 2.0, 2.0], [1.0, 1.0, 1.0], 0.0),     # fit = constant mean
    ([1.0, 2.0, 3.0], [3.0, 2.0, 1.0], [1.0, 1.0, 1.0], -3.0),    # anti-fit → negative (NOT clamped)
    # Zero-weight point must be excluded from BOTH numerator and weighted-mean denominator.
    # w=[1,1,0] with ŷ matching y on the weighted points → R²_w = 1.0 exactly regardless of y_hat[2].
    ([1.0, 2.0, 3.0], [1.0, 2.0, 100.0], [1.0, 1.0, 0.0], 1.0),
    # Partial-weight: ȳ_w = (1·1+1·2+0.5·3)/2.5 = 1.8; SSR = 0; SST = 1·0.64+1·0.04+0.5·1.44 = 1.40 → R² = 1.0
    ([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 1.0, 0.5], 1.0),
])
def test_weighted_r2_pinned(y, y_hat, w, expected_r2_w): ...
```

Negative-R² test is specifically called out because spec forbids clamping; a future "helpful" refactor that clamps would silently break classification (a clamped-to-0 R² would slip a misfit into "Failed" via the low_r2 path instead of flowing through as the true negative value).

**Zero-variance ground truth edge:** `y = [2.0, 2.0, 2.0]` → denominator is 0 → function returns `nan` or a documented sentinel. Test pins that behavior so it is stable.

### 6.7 Classification truth table (R12) — parametrized

Spec table has 5 rows. Test covers all 5 + every boundary.

```python
@pytest.mark.parametrize("r2, n, converged, expected_class, expected_reason", [
    # Row 1: R²≥Good, n in range → Good
    (0.90, 1.0, True, "Good",    None),
    (0.80, 1.0, True, "Good",    None),      # boundary: R² == r2ThresholdGood (inclusive)
    # Row 2: R²≥Good, n out of range → Partial (n_out_of_range downgrade)
    (0.95, 0.4, True, "Partial", None),       # n < nMin → Partial, NOT Failed (spec: "downgrades by one level")
    (0.95, 2.1, True, "Partial", None),
    (0.95, 0.5, True, "Good",    None),       # boundary: n == nMin inclusive
    (0.95, 2.0, True, "Good",    None),       # boundary: n == nMax inclusive
    # Row 3: Failed ≤ R² < Good, n in range → Partial
    (0.50, 1.0, True, "Partial", None),       # boundary: R² == r2ThresholdFailed inclusive
    (0.79, 1.0, True, "Partial", None),
    # Row 4: Failed ≤ R² < Good, n out of range → Failed (n_out_of_range)
    (0.70, 0.3, True, "Failed",  "n_out_of_range"),
    (0.70, 2.5, True, "Failed",  "n_out_of_range"),
    # Row 5: R² < r2ThresholdFailed → Failed (low_r2)
    (0.49, 1.0, True, "Failed",  "low_r2"),
    (-0.5, 1.0, True, "Failed",  "low_r2"),   # negative R² path — classification still correct
    # Did not converge → Failed (convergence_failure), regardless of R²/n
    (0.0,  0.0, False, "Failed", "convergence_failure"),
])
def test_classification_truth_table(r2, n, converged, expected_class, expected_reason): ...
```

**Non-default thresholds:** additional parametrized set runs same table with `r2ThresholdGood=0.95, r2ThresholdFailed=0.7, nMin=0.8, nMax=3.0` (multimeric-antigen scenario). Guards against hardcoded constants.

### 6.8 Failure reason precedence — parametrized

Multiple failure paths can apply simultaneously (e.g., a clonotype below floor AND with non-monotonic signal). Spec implies precedence: insufficient_reads → insufficient_points → non_monotonic_signal → convergence_failure → low_r2/n_out_of_range. Test makes precedence explicit so future refactors don't silently reorder.

```python
@pytest.mark.parametrize("scenario, expected_reason", [
    ("all_below_floor",                    "insufficient_reads"),
    ("some_below_floor_then_under_min",    "insufficient_points"),
    ("enough_points_but_hook",             "non_monotonic_signal"),
    ("hook_but_below_hook_min_reads",      None),                     # proceeds to fit → may be low_r2 or Good
    ("no_hook_but_fit_diverges",           "convergence_failure"),
    ("fit_converges_but_low_r2",           "low_r2"),
    ("fit_converges_good_r2_but_n_oor",    "n_out_of_range"),
    ("good_fit",                           None),
])
def test_failure_reason_precedence(scenario, expected_reason): ...
```

### 6.9 Concentration axis canonicalization (R14) — pinned

```python
# Invariant: floats like 0.001 and 1e-3 may be equal numerically but produce different
# string keys. The pipeline must preserve the user's exact original string.
def test_canonical_axis_roundtrip():
    # Upstream delivers concentrations as strings: ["0", "0.001", "0.01", "0.1", "1"]
    # After fit, output axis keys must be exactly those strings, NOT re-serialized floats.
    ...

def test_canonical_axis_join_does_not_drop_rows():
    # Build reads frame with string concs; run full pipeline; assert every non-c=0 conc
    # appears as an axis key in meanBin output for every fitted clonotype.
    ...

def test_numeric_equal_but_string_different_concs_kept_separate():
    # Pathological: two string concs "1.0" and "1.000" (user-entered differently).
    # Pipeline must NOT merge them — they are distinct axis keys even if numerically equal.
    ...
```

This has historically been a silent failure class in the pframes pipeline. Worth three tests.

### 6.10 Performance guardrail — marked slow

```python
@pytest.mark.slow
@pytest.mark.hot_path
def test_1k_clonotypes_under_wall_clock_budget():
    reads = build_synthetic_reads(n_clonotypes=1_000, n_concs=8, n_bins=4, rng_seed=0)
    t0 = time.perf_counter()
    main.run(reads, params=DEFAULT_PARAMS)
    elapsed = time.perf_counter() - t0
    assert elapsed < 30.0  # generous; alarm for catastrophic regressions only
```

Not a benchmark. Fails only on 10× slowdowns.

---

## 7 — Test data strategy

**Inline (preferred):** most unit tests build a `pl.DataFrame` in-function with 2–6 rows. Each test's dataframe exists to exercise one branch.

**Fixture builders in `conftest.py`:**
- `build_reads(clonotypes, concs, bins, reads_matrix)` — construct canonical long-format reads frame with string-canonical concentrations.
- `build_poisson_reads(true_kd, true_n, n_clonotypes, n_concs, n_bins, seed)` — for regression/sim tests; implements the R-test's spec ("simulate at the read level … not Gaussian noise directly on mean bin").
- `default_params` — frozen dict of R10 defaults for re-use.

**Fixed TSV files in `tests/data/`** (small, for golden-value tests):
- `golden_single_clonotype.tsv` — one clonotype, 8 concentrations including c=0, 4 bins; expected K_D,app baked into test.
- `all_failed.tsv` — every clonotype fails at least one gate; verifies pipeline completes without error.

**What we do NOT do:**
- Do not use production data.
- Do not mock scipy or polars — mock only at the `main.py` argparse boundary (e.g., in `test_io_layer.py`, use `tmp_path` for input files; no monkeypatching of library internals).
- Do not assert on the order in which internal functions are called.

---

## 8 — Open questions to resolve before coding

These are spec ambiguities the plan flags so the implementer doesn't guess:

1. **Hook-effect threshold inclusivity** (§6.3). Spec wording "signal drop at max concentration triggering non_monotonic_signal" doesn't state whether `drop == threshold` flags or not. Plan uses strict `>` (equal doesn't flag). Confirm with spec author.
2. **Weighted R² with zero weighted variance in y** (§6.6). Denominator = 0 edge. Plan returns `nan`; spec is silent. Confirm acceptable — if not, pick a sentinel that `classify` can recognise.
3. **`kdOutOfRange` inclusivity at boundaries**. R14b says "outside [min_concentration, max_concentration]". Plan treats boundary as in-range (closed interval). Confirm.
4. **Failure precedence** (§6.8). Plan pins an order; spec doesn't state one explicitly. Confirm, or add to spec.
5. **Per-clonotype c=0 points that survive the floor but produce `mean_bin = 0`** — do they count toward B? Plan: yes (a valid mean_bin of 0 is a valid baseline observation). Confirm.
6. **Minimum δ value.** Spec references `δ = 0.05` in one place; bounds write `baseline ∈ [0, 0.95]` implying `1 − δ = 0.95`. Plan uses `DELTA = 0.05`. Confirm this is the authoritative value.

---

## 9 — Coverage & review gates

- `uv run pytest --cov=src --cov-report=term-missing tests/` — expect ≥ 90% line coverage on `normalization.py`, `hill_fit.py`, `classify.py`, `floor_filter.py`, `baseline.py`, `hook_effect.py`, `output_build.py`. These are the computational core; gaps here are bugs.
- Lower-priority files (`io_layer.py`, `main.py`) can sit at ≥ 75% — CLI boilerplate and error-path argparse branches are less valuable to test exhaustively.
- Final Phase 4 review (per `python-testing` skill): every test must survive refactoring internals without changing behavior. A test is flagged for rewrite if it asserts on private helpers, call order, or names that don't appear in the spec.

---

## 10 — Out of scope for this plan

- Tengo workflow, model, UI — covered separately.
- Lead Selection integration (R20, R21) — exercised by the block's `test/src/wf.test.ts` end-to-end suite, not by the Python unit tests.
- K_D,app standard error from covariance matrix — explicitly deferred to v2 per README "Out of scope" section; plan does not add tests for it.
- Absolute K_D calibration, multi-experiment comparison, per-lineage stats — all v2+.
