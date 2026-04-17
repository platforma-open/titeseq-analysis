# Test plan — spec bugs from SPEC_REVIEW.md

Behavioral tests that FAIL on `main` today and PASS once each bug is fixed.
All tests use synthetic inline data; no fixtures, no mocks. Run with:

```
cd software/
uv sync && uv run pytest tests/unit/ -v
```

---

## 1. R9b hook effect — missing top-3 half-threshold check

**Bug location:** `software/src/pre_fit.py::detect_hook_effect`
Currently only ranks top-2 and compares `top2 − top1 > threshold`.
Spec requires **both** `top2 − top1 > threshold` **AND** `top3 − top1 > threshold/2`,
with all three top-1/2/3 reads ≥ `hook_effect_min_reads`.

**Target file:** extend `tests/unit/test_pre_fit.py::TestHookEffectBinMode`.

### 1a. Add helper + parametrized cases with 3 concentration points

Replace `_build_fit_points(top2, top1, top2_reads, top1_reads)` with a 3-point
builder and add cases that force the spec's second clause to matter:

```python
def _build_fit_points_3(top3_signal, top2_signal, top1_signal,
                        top3_reads=100, top2_reads=100, top1_reads=100):
    """Three concentrations 1/10/100; signals ordered top-3/top-2/top-1."""
    return pl.DataFrame([
        {"clonotypeKey": "A", "concentrationStr": "1",   "concentration": 1.0,
         "signal": top3_signal, "clonotype_reads_at_conc": top3_reads, "weight": float(top3_reads)},
        {"clonotypeKey": "A", "concentrationStr": "10",  "concentration": 10.0,
         "signal": top2_signal, "clonotype_reads_at_conc": top2_reads, "weight": float(top2_reads)},
        {"clonotypeKey": "A", "concentrationStr": "100", "concentration": 100.0,
         "signal": top1_signal, "clonotype_reads_at_conc": top1_reads, "weight": float(top1_reads)},
    ])


class TestHookEffectBinModeTop3:
    """R9b requires BOTH: top2−top1 > δ AND top3−top1 > δ/2 (δ = 0.2 bin default)."""

    @pytest.mark.parametrize(
        "top3, top2, top1, expected, case",
        [
            # top2-top1 = 0.5 > 0.2 (passes first), top3-top1 = 0.02 ≤ 0.1 (fails second)
            # → spec: NOT flagged (e.g. genuine dose-response that drops only at the very top)
            (2.52, 3.0, 2.5, False, "first_cond_met_second_not_met"),
            # top2-top1 = 0.5 > 0.2 AND top3-top1 = 0.3 > 0.1 → flag
            (2.8, 3.0, 2.5, True, "both_cond_met"),
            # top2-top1 = 0.15 < 0.2 (fails first) → no flag regardless of top-3
            (3.0, 2.95, 2.8, False, "first_cond_not_met"),
            # Edge case: top3 exactly at half-threshold below top1 (0.1 == 0.2/2)
            # spec uses strict >, so should NOT flag
            (2.6, 3.0, 2.5, False, "half_threshold_boundary_equal"),
            # Just above half-threshold on top-3 arm → flag
            (2.601, 3.0, 2.5, True, "half_threshold_boundary_just_over"),
        ],
        ids=lambda x: x if isinstance(x, str) else None,
    )
    def test_bin_mode_top3_clause(self, top3, top2, top1, expected, case):
        df = _build_fit_points_3(top3, top2, top1)
        result = detect_hook_effect(df, bin_mode=True, params=DEFAULT_PARAMS)
        assert bool(result.filter(pl.col("clonotypeKey") == "A")["hook_flag"][0]) is expected, case

    # Spec extends min-reads guard to all three top points.
    @pytest.mark.parametrize(
        "top3_reads, top2_reads, top1_reads, expected",
        [
            (100, 100, 100, True),   # all ≥ 20 → evaluate, drop qualifies
            (10,  100, 100, False),  # top-3 below floor → skip (drop may be noise)
            (100, 10,  100, False),  # top-2 below floor → skip
            (100, 100, 10,  False),  # top-1 below floor → skip
        ],
    )
    def test_min_reads_gate_extends_to_top3(self, top3_reads, top2_reads, top1_reads, expected):
        # Signals chosen so both conditions are satisfied when gate allows
        df = _build_fit_points_3(2.8, 3.0, 2.5,
                                 top3_reads=top3_reads, top2_reads=top2_reads, top1_reads=top1_reads)
        result = detect_hook_effect(df, bin_mode=True, params=DEFAULT_PARAMS)
        assert bool(result.filter(pl.col("clonotypeKey") == "A")["hook_flag"][0]) is expected


class TestHookEffectNoBinModeTop3:
    """Same logic, δ = 0.02."""

    @pytest.mark.parametrize(
        "top3, top2, top1, expected",
        [
            # top2-top1 = 0.03 > 0.02 AND top3-top1 = 0.02 > 0.01 → flag
            (0.04, 0.05, 0.02, True),
            # top2-top1 = 0.03 > 0.02 but top3-top1 = 0.005 < 0.01 → NOT flagged
            (0.025, 0.05, 0.02, False),
        ],
    )
    def test_no_bin_mode_top3_clause(self, top3, top2, top1, expected):
        df = _build_fit_points_3(top3, top2, top1)
        result = detect_hook_effect(df, bin_mode=False, params=DEFAULT_PARAMS)
        assert bool(result.filter(pl.col("clonotypeKey") == "A")["hook_flag"][0]) is expected
```

**Expected result today:** the `first_cond_met_second_not_met` and
`half_threshold_boundary_equal` cases FAIL (current code flags them because it
ignores top-3). The min-reads gate test for `top3_reads=10` also FAILS (current
code only checks top-1/top-2 reads).

---

## 2. R10 DELTA constant — bin mode should use δ = 0.5

**Bug location:** `software/src/constants.py::DELTA` (single value 0.05)
and `software/src/hill_fit.py::fit_one_clonotype` (uses `DELTA` regardless of mode).

**Target file:** extend `tests/unit/test_hill_fit.py::TestHillFitFailureModes`.

```python
class TestDeltaDynamicRangeGate:
    """R10: top − baseline must be ≥ δ_mode. δ_bin = 0.5, δ_no_bin = 0.05."""

    @pytest.mark.parametrize(
        "bin_mode, amplitude, should_converge",
        [
            # Bin mode, δ = 0.5:
            (True,  math.log(0.3), False),  # 0.3 < 0.5 → reject (currently passes — BUG)
            (True,  math.log(0.6), True),   # 0.6 ≥ 0.5 → accept
            (True,  math.log(2.0), True),   # healthy signal
            # No-bin mode, δ = 0.05:
            (False, math.log(0.03), False),  # 0.03 < 0.05 → reject
            (False, math.log(0.08), True),   # 0.08 ≥ 0.05 → accept
        ],
    )
    def test_mode_specific_delta_gate(self, bin_mode, amplitude, should_converge):
        """Generate a clean Hill curve with given amplitude; fit must only
        converge when top − baseline ≥ mode-specific δ."""
        baseline = 1.0 if bin_mode else 0.01
        max_bin_label = 8 if bin_mode else None
        x = LOG_CONCS
        y = hill_truth(x, baseline, amplitude, kd=10.0, n=1.0)
        w = np.ones_like(x) * 100.0

        fit = fit_one_clonotype(x, y, w, baseline_fixed=baseline,
                                bin_mode=bin_mode, max_bin_label=max_bin_label)
        assert fit.converged is should_converge, (
            f"{'bin' if bin_mode else 'no-bin'} amp={math.exp(amplitude):.3f}: "
            f"expected converged={should_converge}, got {fit.converged}"
        )
        if not should_converge:
            assert fit.reason == "convergence_failure"
```

**Expected result today:** the bin-mode `amplitude=log(0.3)` row FAILS — current
code accepts it because `top − baseline = 0.3 ≥ DELTA = 0.05`. After fix
(`DELTA_BIN = 0.5`), it correctly rejects.

---

## 3. R10 amp_lo bounds — solver shouldn't explore below spec δ

**Bug location:** `software/src/hill_fit.py::fit_one_clonotype`
(`amp_lo = math.log(DELTA * 1e-6)` — way below spec bound of `log(δ_mode)`).

This is covered *behaviorally* by test #2: once the post-fit δ gate is
mode-specific, any solver exploration below δ is rejected. A dedicated test
adds a belt-and-braces check:

**Target file:** new test class in `tests/unit/test_hill_fit.py`.

```python
class TestAmplitudeBoundsRespected:
    """Every converged fit must have top − baseline ≥ δ_mode (R10 reparametrization)."""

    @pytest.mark.parametrize("bin_mode, delta_expected", [(True, 0.5), (False, 0.05)])
    def test_converged_top_minus_baseline_ge_delta(self, bin_mode, delta_expected):
        """Well-behaved signal; check the reported top − baseline satisfies spec bound."""
        baseline = 1.0 if bin_mode else 0.01
        max_bin_label = 8 if bin_mode else None
        x = LOG_CONCS
        y = hill_truth(x, baseline, math.log(1.5 if bin_mode else 0.3), kd=10.0, n=1.0)
        w = np.ones_like(x) * 100.0

        fit = fit_one_clonotype(x, y, w, baseline_fixed=baseline,
                                bin_mode=bin_mode, max_bin_label=max_bin_label)
        assert fit.converged is True
        assert (fit.top - fit.baseline) >= delta_expected - 1e-9, (
            f"top − baseline = {fit.top - fit.baseline:.4f} violates δ = {delta_expected}"
        )
```

**Expected result today:** passes on healthy data, but guards the bound for
future refactors. The real catch is test #2.

---

## 4. R5 narrow concentration range warning

**Bug location:** `software/src/io_layer.py::validate_concentration_column`
returns warnings but has no check for `max/min < 10` on non-zero concentrations.

**Target file:** extend `tests/unit/test_io_layer.py::TestConcentrationValidation`.

```python
class TestNarrowConcentrationRangeWarning:
    """R5 guardrail: warn if non-zero concentrations span < 1 order of magnitude."""

    @pytest.mark.parametrize(
        "concs, expect_warning",
        [
            ([1.0, 2.0, 3.0, 5.0, 7.0], True),    # 7/1 = 7 < 10 → warn
            ([1.0, 10.0, 100.0], False),          # 100/1 = 100 ≥ 10 → no warn
            ([1.0, 10.0], False),                 # 10/1 = 10 exactly at boundary → no warn
            ([1.0, 9.99], True),                  # just below → warn
            ([0.0, 1.0, 10.0, 100.0], False),     # 0 excluded from ratio → 100/1 ok
            ([0.0, 1.0, 2.0, 3.0], True),         # 0 excluded → 3/1 < 10 → warn
        ],
    )
    def test_narrow_range_emits_warning(self, concs, expect_warning):
        rows = [
            {"clonotypeKey": "A", "sampleId": f"s{i}",
             "concentrationStr": str(c), "concentration": c, "bin": 1, "reads": 5}
            for i, c in enumerate(concs)
        ]
        df = _mk(rows)
        warnings = validate_concentration_column(df, has_bin=True)
        has_narrow = any("order of magnitude" in w or "narrow" in w.lower() for w in warnings)
        assert has_narrow is expect_warning

    # Edge: a single non-zero concentration — vacuously narrow, should warn.
    def test_single_nonzero_concentration_warns(self):
        df = _mk([
            {"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "0",
             "concentration": 0.0, "bin": 1, "reads": 5},
            {"clonotypeKey": "A", "sampleId": "s2", "concentrationStr": "5",
             "concentration": 5.0, "bin": 1, "reads": 5},
        ])
        warnings = validate_concentration_column(df, has_bin=True)
        assert any("narrow" in w.lower() or "range" in w.lower() for w in warnings)
```

**Expected result today:** all five parametrized cases FAIL (no narrow-range
warning produced). They pass after the validator is extended.

---

## 5. R5 bin × concentration missing-combination warning

**Bug location:** spec section on R5 sub-clause — no implementation anywhere.

Requires a new validator, so the test is paired with a new function.
Proposed signature (for the planner — implementor decides name):

```python
# In io_layer.py:
def validate_bin_concentration_grid(df: pl.DataFrame) -> list[str]:
    """R5: warn if any (bin, concentrationStr) combination present for some
    clonotypes is absent for others."""
```

**Target file:** new `tests/unit/test_io_layer.py::TestBinConcentrationGrid`.

```python
class TestBinConcentrationGrid:
    """R5: warn when (bin, conc) combos are non-uniformly populated across clonotypes."""

    def test_uniform_grid_no_warning(self):
        # Both clonotypes have all (bin, conc) combinations present.
        rows = []
        for clone in ["A", "B"]:
            for conc in ["1", "10"]:
                for b in [1, 2]:
                    rows.append({
                        "clonotypeKey": clone, "sampleId": f"s-{clone}-{conc}-{b}",
                        "concentrationStr": conc, "concentration": float(conc),
                        "bin": b, "reads": 5,
                    })
        warnings = validate_bin_concentration_grid(_mk(rows))
        assert warnings == []

    def test_missing_combo_emits_warning(self):
        # Clonotype A has (bin=1, conc=10); clonotype B lacks that combo.
        rows = [
            {"clonotypeKey": "A", "sampleId": "s1", "concentrationStr": "1",
             "concentration": 1.0, "bin": 1, "reads": 5},
            {"clonotypeKey": "A", "sampleId": "s2", "concentrationStr": "10",
             "concentration": 10.0, "bin": 1, "reads": 5},
            {"clonotypeKey": "B", "sampleId": "s1", "concentrationStr": "1",
             "concentration": 1.0, "bin": 1, "reads": 5},
            # B missing (bin=1, conc=10)
        ]
        warnings = validate_bin_concentration_grid(_mk(rows))
        assert len(warnings) > 0
        assert any("10" in w and "1" in w for w in warnings)
```

---

## Running the new tests

```bash
cd /Users/paulnewling/Desktop/Code/mictx/blocks/titeseq-analysis/software
uv sync
uv run pytest tests/unit/test_pre_fit.py::TestHookEffectBinModeTop3 -v
uv run pytest tests/unit/test_hill_fit.py::TestDeltaDynamicRangeGate -v
uv run pytest tests/unit/test_io_layer.py::TestNarrowConcentrationRangeWarning -v
```

Expected today: §1, §2, §4 FAIL (catching the bugs). §3 passes. §5 cannot
be added until the new validator is implemented.

## Coverage sanity

After adding tests, run:
```bash
uv run pytest --cov=src --cov-report=term-missing tests/unit/
```
and verify `detect_hook_effect`, `fit_one_clonotype`, and
`validate_concentration_column` all reach the new branches.
