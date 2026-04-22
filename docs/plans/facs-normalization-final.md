# FACS Sort-Fraction Normalization — Final Implementation Plan

**Status:** Final — ready to implement
**Block:** `@platforma-open/platforma-open.titeseq-analysis`
**Supersedes:** `facs-normalization.md` (draft), `facs-normalization-final-python-review.md` (folded in)
**Rationale:** `facs-normalization-review.md`
**Spec:** `docs/text/work/projects/affinity-profiling` (external; do not edit)

---

## Summary

Correct the `mean_bin` bias that Tite-Seq's uneven sort yields introduce. The current `compute_mean_bin` weights each FACS bin by its read share only, implicitly assuming every bin caught the same number of cells. In Tite-Seq that is never true — low-antigen concentrations pile cells into the bottom bin, saturation piles them into the top bin — so `K_D,app` drifts with sort yield, not with affinity. The fix is a single multiplicative re-weight by the sort fraction `C_bc / C_c` inside the existing weighted-mean aggregation.

Three decisions shape this plan.

1. **Accept sort fractions as a metadata column**, not a separate CSV upload. One column (`sortFractionColumnRef?: PlRef`) mirrors the existing `binColumnRef` pattern, rides the current metadata-join path through `reads.tsv`, and eliminates the Bin-XYZ string parser, the concentration-index ambiguity, and the `PlFileInput` UI entirely.
2. **Correct `meanBin` in place** and record provenance via a single `pl7.app/titeseq/facsCorrected` annotation. Output schema stays flat — no new column, no dual-signal branch in the Hill fitter.
3. **Keep expression (c-Myc) normalization explicitly out of scope.** The field does not divide by expression (Adams, Mora, Walczak, Kinney 2016); the Hill amplitude absorbs display differences. Call this out so readers do not mistake absence for oversight.

With those three decisions baked in, the implementation touches five packages (software, workflow, model, UI, block), adds one kwarg to one Python function, one anchor to the Tengo bundle, one field to `BlockArgs`, and one dropdown to the settings drawer. The rest is tests and release hygiene.

---

## 1. The Bias The Fix Removes

The current `compute_mean_bin` in `software/src/normalization.py` computes:

```
freq_{s,b,c} = reads_{s,b,c} / R_{b,c}          # read share of variant s in bin b at conc c
mean_bin(s,c) = Σ_b b · freq / Σ_b freq
```

This is valid only when every bin receives the same cell count at every concentration. Tite-Seq inverts that assumption by construction. The correct formula follows from Adams, Mora, Walczak, Kinney (eLife 2016; e23156), Appendix 5, eq. A3:

```
r_{s,b,c} ≈ (R_bc / C_bc) · C_c · P_s · p(b | s, c)
```

Rearranging and normalizing over bins yields the probability that a variant-`s` cell at concentration `c` lands in bin `b`:

```
p(b | s, c) = [freq_{s,b,c} · (C_bc / C_c)] / Σ_b [freq_{s,b,c} · (C_bc / C_c)]
```

where `C_bc / C_c` is the sort fraction — the user's wet-lab-measured cell yield per bin per concentration. The corrected signal is `mean_bin(s,c) = Σ_b b · p(b | s, c)`. It reduces to today's formula exactly when every bin catches an equal cell fraction, and diverges from it otherwise. The divergence is not small: for a typical 4-bin Tite-Seq sort at saturating concentration, ~80% of cells land in bin 4, so the uncorrected formula over-weights bin 4 by roughly 3×.

---

## 2. Scope

**In scope.** Accept a `sort_fraction` per-sample metadata column. Re-weight per-bin read frequencies by `sort_fraction` before the weighted mean. Fall back to today's formula bit-exactly when the column is absent. Annotate the output to record which path ran.

**Out of scope, named deliberately.**

- **Expression (c-Myc) normalization** — §9a.
- **Bin-MFI log-fluorescence reconstruction** — §9b. Canonical Tite-Seq uses `log μ_b` bin centers rather than integer labels; integer labels are the block's documented simplification.
- **Poisson MLE joint Kd/Hill fit** (Adams eq. A4) — requires gate boundaries the user has not provided.
- **Dirichlet-multinomial Bayesian error model** (Phillips 2021) — principled future upgrade for tight CIs on low-count variants.

---

## 3. Input Data Channel

The user adds one column to the metadata TSV. No file upload, no new code path.

**User-facing format.** One row per sample, value = `C_bc / C_c` for that sample's `(concentration, bin)` pair. Values sum to 1.0 per concentration.

| sample_id | concentration_M | bin | sort_fraction |
|---|---|---|---|
| S001 | 1e-06 | 1 | 0.5068 |
| S002 | 1e-06 | 2 | 0.1176 |
| S003 | 1e-06 | 3 | 0.1492 |
| S004 | 1e-06 | 4 | 0.2264 |

**Migration from the current `facs_df.csv` format.** A one-time reshape (ship this snippet in the PR description):

```python
import polars as pl

facs = pl.read_csv("facs_df.csv")  # columns: Bin, Fraction
facs = facs.with_columns([
    pl.col("Bin").str.extract(r"Bin-(\d{2})", 1).cast(pl.Int64).alias("conc_point"),
    pl.col("Bin").str.extract(r"Bin-\d{2}(\d)", 1).cast(pl.Int64).alias("bin"),
]).rename({"Fraction": "sort_fraction"})

metadata = pl.read_csv("titeseq_metadata.tsv", separator="\t")  # has conc_point, bin
out = metadata.join(facs, on=["conc_point", "bin"], how="left")
out.write_csv("titeseq_metadata_with_facs.tsv", separator="\t")
```

**Validation invariants (Python, `io_layer.py`).**

- All values in `[0, 1]`.
- Per concentration, `|Σ_b sort_fraction − 1| < 1e-3`.
- Every `(sample_id)` with reads has a non-null `sort_fraction`.
- No silent renormalization. Violations raise `InputValidationError` with a concrete message naming the offending concentration.

---

## 4. Data Flow

```
┌──────────────────────────────────────────────────────────────┐
│  User: adds `sort_fraction` column to titeseq_metadata.tsv   │
│  → imports via their usual metadata-ingestion block          │
└────────────┬─────────────────────────────────────────────────┘
             │  PColumn [sampleId] → Double
             ▼
┌──────────────────────────────────────────────────────────────┐
│  Model (model/src/index.ts)                                   │
│    sortFractionColumnRef?: PlRef                              │
│    output: sortFractionColumnOptions                          │
│    validation: requires binColumnRef                          │
└────────────┬─────────────────────────────────────────────────┘
             │
             ▼
┌──────────────────────────────────────────────────────────────┐
│  UI (SettingsDrawer.vue)                                      │
│    PlDropdownRef, gated on binColumnRef present               │
└────────────┬─────────────────────────────────────────────────┘
             │
             ▼
┌──────────────────────────────────────────────────────────────┐
│  Workflow (main.tpl.tengo)                                    │
│    bb.addAnchor("sortFraction", args.sortFractionColumnRef)   │
│    join onto reads.tsv as column "sort_fraction"              │
│    exec arg: --sort-fraction-column sort_fraction             │
│    annotation on meanBin: pl7.app/titeseq/facsCorrected       │
└────────────┬─────────────────────────────────────────────────┘
             │
             ▼
┌──────────────────────────────────────────────────────────────┐
│  Python (software/src/)                                       │
│    main.py:          --sort-fraction-column → run(...)        │
│    io_layer.py:      validate_sort_fraction                   │
│    normalization.py: compute_mean_bin(..., sort_fraction_col) │
│    pipeline.py:      log branch; unchanged return shape       │
└──────────────────────────────────────────────────────────────┘
```

---

## 5. Python Software Changes

### 5a. `software/src/constants.py`

Add the constant with an explicit cross-language contract comment. The Tengo workflow emits the same string into `reads.tsv` and passes it via `--sort-fraction-column`; the two must not drift apart.

```python
# Column name emitted by the Tengo workflow in reads.tsv when
# sortFractionColumnRef is bound. Must match the string passed to
# --sort-fraction-column. Do not change without a simultaneous change
# to workflow/src/main.tpl.tengo.
COL_SORT_FRACTION = "sort_fraction"
```

### 5b. `software/src/normalization.py` — extend `compute_mean_bin`

One optional keyword-only kwarg. Legacy path is a special case (`sort_fraction_col=None`) and must stay bit-exact — that is the regression gate. Keyword-only prevents a future caller from silently disabling the correction by a positional slip.

```python
def compute_mean_bin(
    reads: pl.DataFrame,
    *,
    sort_fraction_col: str | None = None,
) -> pl.DataFrame:
    """Mean bin per (clonotype, concentration).

    Two modes:
    * sort_fraction_col=None (legacy): mean_bin = Σ_b b·freq / Σ_b freq.
    * sort_fraction_col="<col>": Adams, Mora, Walczak, Kinney 2016 eq. A3
      correction; mean_bin = Σ_b b·freq·frac / Σ_b freq·frac. Column must
      be present in `reads` and pre-validated by validate_sort_fraction.
    """
    with_freq = (
        reads.with_columns(
            pl.col(COL_READS).sum().over([COL_CONC_STR, COL_BIN]).alias("depth"),
        )
        .with_columns(
            pl.when(pl.col("depth") > 0)
              .then(pl.col(COL_READS) / pl.col("depth"))
              .otherwise(0.0)
              .alias(FREQ)
        )
    )

    weight = (
        pl.col(FREQ) * pl.col(sort_fraction_col)
        if sort_fraction_col is not None
        else pl.col(FREQ)
    )

    return (
        with_freq.with_columns(weight.alias("w"))
        .group_by([COL_CLONOTYPE, COL_CONC_STR, COL_CONC_VAL])
        .agg(
            (pl.col(COL_BIN).cast(pl.Float64) * pl.col("w")).sum().alias("num"),
            pl.col("w").sum().alias("den"),
            pl.col(COL_READS).sum().alias(CLONOTYPE_READS_AT_CONC),
        )
        .with_columns(
            pl.when(pl.col("den") > 0)
              .then(pl.col("num") / pl.col("den"))
              .otherwise(None)
              .alias(MEAN_BIN)
        )
        .drop(["num", "den"])
    )
```

Extend `normalize` with a keyword-only kwarg for the same reason:

```python
def normalize(
    reads: pl.DataFrame,
    bin_mode: bool,
    *,
    sort_fraction_col: str | None = None,
) -> pl.DataFrame:
    if bin_mode:
        sig = compute_mean_bin(reads, sort_fraction_col=sort_fraction_col)
        return sig.rename({MEAN_BIN: SIGNAL})
    return compute_frequency_signal(reads)
```

### 5c. `software/src/io_layer.py` — add validator

```python
def validate_sort_fraction(reads: pl.DataFrame, col: str) -> None:
    if col not in reads.columns:
        raise InputValidationError(f"Sort-fraction column '{col}' missing from reads")

    vals = reads[col]
    if vals.is_null().any():
        raise InputValidationError(f"Sort-fraction column '{col}' has null values")
    if (vals < 0).any() or (vals > 1).any():
        raise InputValidationError(f"Sort-fraction values must be in [0, 1]; column '{col}'")

    sums = reads.group_by(COL_CONC_STR).agg(pl.col(col).sum().alias("s"))
    bad = sums.filter((pl.col("s") - 1.0).abs() > 1e-3)
    if bad.height > 0:
        violations = bad.select([COL_CONC_STR, "s"]).to_dicts()
        raise InputValidationError(
            f"sort_fraction values must sum to 1.0 per concentration "
            f"(tolerance 1e-3). Violations: {violations}"
        )
```

Only the violating rows appear in the message — not the full per-concentration dump — so the user sees exactly which concentration to fix.

Call once in `pipeline.run` when `sort_fraction_col` is set, after the reads frame is loaded.

### 5d. `software/src/pipeline.py` — thread kwarg, log branch

```python
def run(
    reads: pl.DataFrame,
    *,
    params: FitParams,
    target_antigen: str | None,
    antigen_column_ref: str | None,
    sort_fraction_col: str | None = None,
) -> PipelineOutputs:
    if sort_fraction_col is not None:
        validate_sort_fraction(reads, sort_fraction_col)
        n_conc = reads.select(COL_CONC_STR).n_unique()
        log.info(
            f"sort_fraction validated: {n_conc} concentrations, "
            f"all sums within 1e-3 of 1.0"
        )
        log.info(f"Mean-bin correction: FACS-weighted (column '{sort_fraction_col}')")
    else:
        log.info("Mean-bin correction: uncorrected")

    # ... existing flow, thread sort_fraction_col into normalize(...)
```

Logging the validator pass explicitly prevents silent validator regressions from slipping past QA review.

### 5e. `software/src/main.py` — CLI

```python
parser.add_argument(
    "--sort-fraction-column",
    default=None,
    help="Column name in reads.tsv carrying sort fraction C_bc/C_c per sample. "
         "Absent ⇒ no FACS correction.",
)
# ...
outputs = run(reads, ..., sort_fraction_col=args.sort_fraction_column)
```

### 5f. `software/src/output_build.py`

No schema change. `meanBin` column already exists; the annotation comes from the workflow layer.

### 5g. `software/package.json`

Bump `version`. CI republishes the Python binary only on version change.

---

## 6. Workflow (Tengo) Changes

### 6a. `workflow/src/main.tpl.tengo` — bundle builder

Add a fifth anchor in `wf.prepare`:

```go
wf.prepare(func(args) {
	bb := wf.createPBundleBuilder()
	bb.ignoreMissingDomains()
	bb.addAnchor("abundance", args.abundanceRef)
	bb.addAnchor("concentration", args.concentrationColumnRef)

	if !is_undefined(args.binColumnRef) {
		bb.addAnchor("bin", args.binColumnRef)
	}
	if !is_undefined(args.antigenColumnRef) {
		bb.addAnchor("antigen", args.antigenColumnRef)
	}
	if !is_undefined(args.sortFractionColumnRef) {
		bb.addAnchor("sortFraction", args.sortFractionColumnRef)
	}

	return { columns: bb.build() }
})
```

### 6b. `wf.body` — join column into reads TSV

Mirror the existing `bin` column join in the `tsvFileBuilder`:

```go
hasSortFraction := !is_undefined(args.sortFractionColumnRef)
if hasSortFraction {
	sortFractionData := columns.getData("sortFraction")
	tsvBuilder.addColumn("sort_fraction", sortFractionData)
}
```

(The precise call is whatever idiom `tsvFileBuilder` uses for `bin`. Copy that pattern; do not invent a new one.)

### 6c. `exec.builder` — pass the argument

```go
execBuilder := exec.builder().
	software(fitCurvesSw).
	mem("16GiB").
	cpu(4).
	writeFile("params.json", json.encode(paramsPayload)).
	addFile("reads.tsv", readsTsv).
	arg("--reads").arg("reads.tsv").
	arg("--params").arg("params.json").
	arg("--out-per-clonotype").arg("per_clonotype.tsv").
	arg("--out-mean-bin").arg("mean_bin.tsv").
	arg("--out-fitted-mean-bin").arg("fitted_mean_bin.tsv").
	saveFile("per_clonotype.tsv").
	saveFile("mean_bin.tsv").
	saveFile("fitted_mean_bin.tsv").
	saveStdoutStream()

if hasSortFraction {
	execBuilder = execBuilder.arg("--sort-fraction-column").arg("sort_fraction")
}
```

### 6d. Annotate `meanBin` with provenance

When constructing the `meanBin` PColumn spec:

```go
facsCorrected := "false"
if hasSortFraction {
	facsCorrected = "true"
}

meanBinAnnotations := {
	"pl7.app/label": meanBinLabel,
	"pl7.app/table/orderPriority": "90000",
	"pl7.app/trace": trace.toJson(),
	"pl7.app/titeseq/facsCorrected": facsCorrected
}
```

Downstream Hill fit and Graph Maker see the same column shape; the annotation propagates through `pSpec.makeTrace`.

### 6e. `paramsPayload` — audit log

```go
paramsPayload.facs_correction_active = hasSortFraction
```

Appears in the Python log and in block state via `get_block_state` transforms.

### 6f. Version pin

Bump the version in the `importSoftware` line to match `software/package.json`:

```go
fitCurvesSw := assets.importSoftware("@platforma-open/platforma-open.titeseq-analysis.software:fit-curves")
```

The catalog entry (not shown here) carries the pinned version — update it.

---

## 7. Model (TypeScript) Changes

### 7a. `model/src/index.ts` — add the field

```typescript
export type BlockArgs = {
  abundanceRef?: PlRef;
  concentrationColumnRef?: PlRef;
  binColumnRef?: PlRef;
  antigenColumnRef?: PlRef;
  sortFractionColumnRef?: PlRef;   // NEW
  targetAntigen?: string;
  // ... remaining fields unchanged
};
```

`DataModelBuilder.init()`:

```typescript
.init(() => ({
  abundanceRef: undefined,
  concentrationColumnRef: undefined,
  binColumnRef: undefined,
  antigenColumnRef: undefined,
  sortFractionColumnRef: undefined,   // NEW
  targetAntigen: undefined,
  // ...
}))
```

### 7b. `.args<BlockArgs>` — one validation

```typescript
if (data.sortFractionColumnRef !== undefined && data.binColumnRef === undefined) {
  throw new Error("Sort-fraction column requires a bin column");
}
```

Value-level validation (range, sum-to-one) lives in the Python software — enforcing it in TS would duplicate logic and require materializing the column.

### 7c. Options output

```typescript
.output("sortFractionColumnOptions", (ctx) =>
  ctx.resultPool.getOptions([{
    axes: [{ name: "pl7.app/sampleId" }],
    annotations: {},
  }])
)
```

Filter to numerical columns in the UI layer if `getOptions` does not already.

### 7d. Provenance read-through

```typescript
.output("facsCorrectionActive", (ctx) => {
  const spec = ctx.outputs?.resolve("meanBin")?.getPColumnSpec?.();
  return spec?.annotations?.["pl7.app/titeseq/facsCorrected"] === "true";
})
```

Surfaces to the UI badge.

---

## 8. UI (Vue) Changes

### 8a. `ui/src/components/SettingsDrawer.vue`

Under Inputs, immediately after the bin column dropdown:

```vue
<PlDropdownRef
  v-if="app.model.data.binColumnRef !== undefined"
  v-model="app.model.data.sortFractionColumnRef"
  :options="app.model.outputs.sortFractionColumnOptions"
  label="FACS sort fraction (optional)"
  clearable
>
  <template #tooltip>
    Optional. A per-sample numerical metadata column whose value is the
    fraction of cells sorted into that sample's (concentration, bin) —
    C_bc / C_c in Adams, Mora, Walczak, Kinney 2016. Values must sum to 1
    per concentration. When supplied, Mean Bin is corrected for FACS sort
    yield.
  </template>
</PlDropdownRef>
```

### 8b. Header badge

In `ui/src/components/PageHeader.vue`:

```vue
<span v-if="app.model.outputs.facsCorrectionActive" class="badge badge-success">
  FACS-corrected
</span>
```

---

## 9. Explicit Deferrals

### 9a. Expression (c-Myc) normalization

Adams, Mora, Walczak, Kinney 2016 does not divide `mean_bin` by the expression channel. The c-Myc signal is used for filtering non-displayers and for QC. The Hill amplitude `A_v` absorbs per-variant display differences because `K_D,app` is a half-saturation parameter invariant to amplitude.

The spec's R6 limitation note documents the residual bias this leaves in libraries with large display variation (CDR walks with hydrophobic mutations). The right future fix is an expression-channel filter, not a per-variant divisor:

- New metadata column `expressionColumnRef?: PlRef` (per-variant mean display signal).
- Python pre-fit filter: drop variants with mean expression below a threshold.
- Same PColumn-anchor pattern as `sortFractionColumnRef`. No fitter changes.

### 9b. Bin-MFI log-fluorescence signal

Canonical Tite-Seq uses geometric-mean log-fluorescence bin centers `μ_b` rather than integer labels. Integer labels are the block's documented simplification (`pcolumn-spec.md`). Upgrading swaps `b` for `log μ_b` in the same weighted-mean formula, so sort-fraction correction composes cleanly. Requires a second metadata channel (bin-center fluorescence per concentration). Separate ticket.

### 9c. Poisson MLE joint Kd/Hill fit (Adams eq. A4)

Requires per-concentration per-bin gate boundaries the user has not provided. Not blocked by this plan.

### 9d. Dirichlet-multinomial Bayesian model (Phillips 2021)

Principled tight-CI upgrade for low-count variants. Defer until users ask.

---

## 10. Output Schema Delta

No new columns. One new annotation on the existing `meanBin` PColumn:

| Annotation | Values | Meaning |
|---|---|---|
| `pl7.app/titeseq/facsCorrected` | `"true"` \| `"false"` | Whether sort-fraction correction ran |

Spec R14 (`meanBin` PColumn with axes `[clonotypeKey][concentration]`) is satisfied unchanged. The annotation propagates via `pSpec.makeTrace` like every other PColumn annotation.

---

## 11. Edge Cases

| Case | Behaviour |
|---|---|
| `sortFractionColumnRef` undefined | Legacy path. No annotation change from pre-feature output (annotation value = `"false"`). |
| `sort_fraction = 0` for a `(conc, bin)` | Bin contributes nothing to numerator or denominator. If every bin at a concentration is zero, `mean_bin` is null there and the existing `minConcentrationPoints` filter applies. |
| Sort fractions sum to 0.98 or 1.02 at some concentration | Rejected at ingestion (`InputValidationError`). No silent renormalization. |
| `sort_fraction` null for a sample that has reads | Rejected at ingestion. The column must be complete if supplied. |
| `sortFractionColumnRef` set, `binColumnRef` unset | TS model validation error — correction applies only in bin mode. |
| Negative or `>1` sort fraction | Rejected at ingestion. |
| User supplies the column but not all samples share its domain | `getOptions` filters the dropdown to valid columns; the anchor resolver surfaces a clear error if ambiguity remains. |

---

## 12. Testing

All tests assert on return values, shapes, raised exceptions, or pinned outputs. None mock internal functions, assert call counts, or touch private attributes.

### 12a. Unit — `software/tests/unit/test_normalization.py`

**Shared fixture.** Extend the existing `_build_single_clonotype` helper with one optional parameter. Default `None` leaves legacy tests untouched; FACS tests reuse the same construction path.

```python
def _build_single_clonotype(
    reads_per_bin: list[int],
    depth_per_bin: list[int],
    sort_fraction_per_bin: list[float] | None = None,
) -> pl.DataFrame:
    ...
    if sort_fraction_per_bin is not None:
        row["sort_fraction"] = sort_fraction_per_bin[j]
    ...
```

**Test cases.**

- `test_compute_mean_bin_legacy_bit_exact` — `compute_mean_bin(reads, sort_fraction_col=None)` equals a frozen reference CSV byte-for-byte. **Regression gate.**
- `test_compute_mean_bin_uniform_fractions` — `sort_fraction = 1/|bins|` uniformly ⇒ output equals legacy within float tolerance.
- `test_compute_mean_bin_skewed_fractions_pinned` — exact numeric values, not direction. Parametrize over several skews. One anchored case:

  ```python
  # 2 bins, freq=[0.1, 0.1], fraction=[0.8, 0.2].
  # Weighted numerator = 1·(0.1·0.8) + 2·(0.1·0.2) = 0.12.
  # Weighted denominator = 0.1·0.8 + 0.1·0.2 = 0.10.
  # Corrected mean_bin = 1.2. Uncorrected would be 1.5.
  def test_compute_mean_bin_skewed_fractions_pinned(self):
      reads = _build_single_clonotype(
          reads_per_bin=[10, 10],
          depth_per_bin=[100, 100],
          sort_fraction_per_bin=[0.8, 0.2],
      )
      out = compute_mean_bin(reads, sort_fraction_col="sort_fraction")
      assert out["mean_bin"][0] == pytest.approx(1.2)
  ```

- `test_compute_mean_bin_zero_fraction_bin` — bin with `sort_fraction = 0` excluded from both numerator and denominator, including when it has nonzero reads.
- `test_compute_mean_bin_per_concentration_independence` — two concentrations, same clonotype, different sort skews per concentration. Guards against the correction leaking across the `group_by([clonotype, concentration])` key. Single-concentration tests cannot catch this class of bug.
- `test_normalize_dispatch_legacy` — `normalize(reads, bin_mode=True, sort_fraction_col=None)` bit-exact vs. today.

### 12b. Unit — `software/tests/unit/test_io_layer.py`

- Happy path — in-range, sum ≈ 1 per concentration.
- Out-of-range — negative values, values `> 1`.
- Sum-violation — fractions summing to 0.9 or 1.1 at one concentration; assert the error message names the offending concentration:

  ```python
  def test_validate_sort_fraction_names_offending_concentration(self):
      # User needs to know WHICH concentration is misconfigured.
      df = _mk_reads_with_sort_fraction(
          concentrations=["1", "10"],
          sort_fractions_by_conc={
              "1": [0.25, 0.25, 0.25, 0.25],
              "10": [0.5, 0.3, 0.1, 0.0],  # sums to 0.9
          },
      )
      with pytest.raises(InputValidationError, match="10"):
          validate_sort_fraction(df, "sort_fraction")
  ```

- Missing-column — validator raises before `compute_mean_bin` ever sees the frame.
- Null-propagation — a join miss produces `None` in the column. Pin the expected behavior: validator rejects, reads are never silently lost:

  ```python
  def test_validate_sort_fraction_rejects_null(self):
      # A join miss produces null; the user should see it at validation,
      # not silently lose reads through polars' null-skipping aggregations.
      df = ...  # row with sort_fraction = None
      with pytest.raises(InputValidationError, match="null"):
          validate_sort_fraction(df, "sort_fraction")
  ```

### 12c. Integration — `software/tests/integration/test_bin_mode_pipeline.py`

Add `test_pipeline_with_sort_fraction`. Use a BOSI-shaped fixture (13 concentrations × 4 bins × ~5 clonotypes). Assert:

- `mean_bin.tsv` shape matches uncorrected run.
- Hill fit succeeds.
- `K_D,app` differs from uncorrected run in a direction consistent with the fixture's injected skew.

### 12d. Regression — `software/tests/regression/test_synthetic_titeseq.py`

Extend the synthetic generator to simulate skewed sort yields. Pin thresholds — a vague "degrades measurably" passes silently when the degradation is small.

```python
# Synthetic ground-truth K_D injected per clonotype. Success =
# recovered K_D within 10% of injected. Corrected path must clear 90%;
# uncorrected path must fall below 70%. A gap < 20 points means the
# fixture skew is too mild to exercise the bias — strengthen it.
def test_synthetic_recovery_corrected_vs_uncorrected(self):
    corrected_rate = run_and_measure(apply_facs=True)
    uncorrected_rate = run_and_measure(apply_facs=False)
    assert corrected_rate >= 0.90
    assert uncorrected_rate < 0.70
    assert corrected_rate - uncorrected_rate > 0.20
```

### 12e. Block integration — `test/src/wf.test.ts`

End-to-end test with `sortFractionColumnRef` bound to a metadata column in the test fixtures. Assert:

- Block runs to completion.
- `meanBin` PColumn present.
- `pl7.app/titeseq/facsCorrected` annotation = `"true"`.
- `facsCorrectionActive` output = `true`.

### 12f. Coverage — expected branch deltas

Run `uv run pytest --cov=src --cov-branch --cov-report=term-missing tests/` after implementation. Target: no uncovered behavioral branches in `normalization.py` or `io_layer.py`.

| Module / function | New branches | Covered by |
|---|---|---|
| `normalization.py::compute_mean_bin` | 2 (ternary on `sort_fraction_col`) | legacy bit-exact + skewed-fractions pinned |
| `normalization.py::normalize` | 1 (kwarg threading) | `test_normalize_dispatch_legacy` |
| `io_layer.py::validate_sort_fraction` | 5 (happy, missing col, out-of-range, sum violation, null) | §12b cases |
| `pipeline.py::run` | 1 (log branch) | Coverage-exempt (logging-only) |
| `main.py::parser` | 1 (CLI arg) | Integration-level; unit coverage not expected |

---

## 13. Changeset

`.changeset/facs-sort-fraction-correction.md`:

```markdown
---
'@platforma-open/platforma-open.titeseq-analysis.software': minor
'@platforma-open/platforma-open.titeseq-analysis.workflow': minor
'@platforma-open/platforma-open.titeseq-analysis.model': minor
'@platforma-open/platforma-open.titeseq-analysis.ui': minor
'@platforma-open/platforma-open.titeseq-analysis': minor
---

Add optional FACS sort-fraction correction to Mean Bin.

Users can now supply a per-sample sort_fraction metadata column
(C_bc/C_c from Adams, Mora, Walczak, Kinney 2016). When present, Mean Bin
is sort-yield-corrected and the output carries the annotation
pl7.app/titeseq/facsCorrected="true". Absent the column, behaviour is
bit-exact with the prior release.
```

All five packages bump at `minor` — minor changes across every layer require an explicit minor bump on the root block package (Changesets only propagates patch automatically).

---

## 14. Implementation Checklist

**Preparation**
- [ ] Verify `@platforma-sdk/block-tools` in `pnpm-workspace.yaml` matches latest npm (`npm view @platforma-sdk/block-tools version`). Update if behind.
- [ ] Confirm clean working tree (`git status`).
- [ ] Run `git fetch && git status` — do not work on stale branches.

**Python software**
- [ ] `src/constants.py` — `COL_SORT_FRACTION`.
- [ ] `src/normalization.py` — `sort_fraction_col` kwarg on `compute_mean_bin` and `normalize`.
- [ ] `src/io_layer.py` — `validate_sort_fraction`.
- [ ] `src/pipeline.py` — thread kwarg; log branch.
- [ ] `src/main.py` — `--sort-fraction-column`.
- [ ] `package.json` — bump version.
- [ ] `tests/unit/test_normalization.py` — 5 new cases (§12a).
- [ ] `tests/unit/test_io_layer.py` — 5 new cases (§12b).
- [ ] `tests/integration/test_bin_mode_pipeline.py` — `test_pipeline_with_sort_fraction`.
- [ ] `tests/regression/test_synthetic_titeseq.py` — skewed-yield path.

**Workflow**
- [ ] `src/main.tpl.tengo` — `sortFraction` anchor.
- [ ] `src/main.tpl.tengo` — join `sort_fraction` into `reads.tsv`.
- [ ] `src/main.tpl.tengo` — pass `--sort-fraction-column sort_fraction` when anchor present.
- [ ] `src/main.tpl.tengo` — `pl7.app/titeseq/facsCorrected` annotation on `meanBin` spec.
- [ ] `src/main.tpl.tengo` — `paramsPayload.facs_correction_active` for audit.
- [ ] Bump `importSoftware` version pin to match software `package.json`.

**Model**
- [ ] `src/index.ts` — `sortFractionColumnRef?: PlRef` in `BlockArgs`.
- [ ] `src/index.ts` — initialize to `undefined` in `DataModelBuilder.init`.
- [ ] `src/index.ts` — validation: sort-fraction requires bin.
- [ ] `src/index.ts` — output `sortFractionColumnOptions`.
- [ ] `src/index.ts` — output `facsCorrectionActive` from annotation.

**UI**
- [ ] `components/SettingsDrawer.vue` — `PlDropdownRef` gated on bin column.
- [ ] `components/PageHeader.vue` — FACS-corrected badge.

**Block integration**
- [ ] `test/src/wf.test.ts` — end-to-end with sort-fraction column.
- [ ] `test/fixtures/` — metadata fixture with `sort_fraction` column.

**Python test verification**
- [ ] `cd blocks/titeseq-analysis/software && uv sync`.
- [ ] `uv run pytest tests/ -v` — all pass.
- [ ] `uv run pytest --cov=src --cov-branch --cov-report=term-missing tests/` — no uncovered behavioral branches in `normalization.py` or `io_layer.py`.

**Release**
- [ ] `.changeset/facs-sort-fraction-correction.md` — per §13.
- [ ] `pnpm install` — lock file updated.
- [ ] `pnpm run build:dev` — local build (user preference).
- [ ] MCP dev cycle in Platforma TiteSeq project: `update_block` → `run_block` → `await_block_done` → `get_block_state` with transform on `meanBin` and `K_D,app`.
- [ ] `git status` clean; `pnpm-lock.yaml` included if `pnpm-workspace.yaml` touched.
- [ ] `git diff --name-only origin/main..HEAD` confirms exactly the intended files.
- [ ] Code comment in `normalization.py` pointing at this file (spec-deviation convention).

---

## 15. References

**Primary source (verified):**

1. Adams RM, Mora T, Walczak AM, Kinney JB. "Measuring the sequence-affinity landscape of antibodies with massively parallel titration curves." *eLife* 5:e23156 (2016). https://elifesciences.org/articles/23156 — Appendix 5, eq. A3 (sort-fraction correction derivation).

**Supporting:**

2. Starr TN et al. "Deep mutational scanning of SARS-CoV-2 receptor binding domain reveals constraints on folding and ACE2 binding." *Cell* 182:1295 (2020). https://doi.org/10.1016/j.cell.2020.08.012 — production-scale 4-bin titration.
3. Phillips AM et al. "Hierarchical multi-replicate Bayesian analysis of high-throughput titration curves." *eLife* 10:e71393 (2021) — Bayesian error model (future milestone).
4. Bloom lab `dms_variants`. https://jbloomlab.github.io/dms_variants/ — Python reference implementation (uses bin-MFI rather than integer labels).

**Internal:**

5. `docs/text/work/projects/affinity-profiling/pcolumn-spec.md` — integer bin labels, R6 global baseline, R14 output schema.
6. `docs/text/work/projects/affinity-profiling/README.md` — scope ("metadata-only inputs").
7. `blocks/titeseq-analysis/docs/plans/facs-normalization.md` — prior draft (superseded).
8. `blocks/titeseq-analysis/docs/plans/facs-normalization-review.md` — review decisions driving this plan.
9. `blocks/titeseq-analysis/software/src/normalization.py:compute_mean_bin` — current unweighted implementation.
10. `blocks/titeseq-analysis/workflow/src/main.tpl.tengo` — bundle-builder anchor pattern; PColumn spec emission.
