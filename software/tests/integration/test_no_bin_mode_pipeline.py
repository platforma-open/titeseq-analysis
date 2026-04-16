"""End-to-end no-bin mode pipeline tests (R7b signal = frequency)."""

from __future__ import annotations

import math

import polars as pl
import pytest

from constants import DEFAULT_PARAMS
from pipeline import run


def _build_no_bin_reads_for_hill(
    clonotype_reads_at_conc: dict[str, list[int]],
    total_reads_at_conc: list[int],
    conc_strs: list[str],
    concs: list[float],
) -> pl.DataFrame:
    """For each concentration, each clonotype k has a reads count; the difference
    between sum(clonotype reads) and total_reads_at_conc is added under a filler clonotype."""
    rows = []
    for i, (cs, c) in enumerate(zip(conc_strs, concs)):
        total = total_reads_at_conc[i]
        used = 0
        for k, reads_list in clonotype_reads_at_conc.items():
            r = reads_list[i]
            used += r
            rows.append(
                {
                    "clonotypeKey": k,
                    "sampleId": f"s_c{i}",
                    "concentrationStr": cs,
                    "concentration": c,
                    "reads": r,
                }
            )
        filler = total - used
        if filler > 0:
            rows.append(
                {
                    "clonotypeKey": "Filler",
                    "sampleId": f"s_c{i}",
                    "concentrationStr": cs,
                    "concentration": c,
                    "reads": filler,
                }
            )
    return pl.DataFrame(rows)


class TestNoBinModePipeline:
    # No-bin mode: signal = freq; fit a known Hill response.
    def test_good_clonotype_recovered_no_bin(self):
        concs = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0]
        conc_strs = [str(c) for c in concs]
        baseline = 0.001
        amplitude = math.log(0.3)
        true_kd = 10.0
        true_n = 1.0

        def hill(c):
            top = baseline + math.exp(amplitude)
            return baseline + (top - baseline) * (c**true_n) / (true_kd**true_n + c**true_n)

        total = 100_000
        # Each concentration: clonotype "G" takes freq*total reads (rounded), filler gets rest.
        reads_g = [round(hill(c) * total) for c in concs]
        df = _build_no_bin_reads_for_hill(
            {"G": reads_g},
            total_reads_at_conc=[total] * len(concs),
            conc_strs=conc_strs,
            concs=concs,
        )
        out = run(df, params=DEFAULT_PARAMS)
        pc = out["per_clonotype"].filter(pl.col("clonotypeKey") == "G")
        assert pc["affinityClass"][0] in {"Good", "Partial"}
        assert pc["kd"][0] == pytest.approx(10.0, rel=0.3)
