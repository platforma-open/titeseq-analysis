"""Deterministic synthetic MiXCR-format TSV fixture generator for block tests.

Emits per-sample MiXCR-format TSVs (headers: readCount, nSeqCDR3, aaSeqCDR3,
bestVGene, bestJGene) plus a manifest.json describing samples, metadata, and
expected per-clonotype outcomes. Block tests import these via the S&D +
import-vdj chain (format='mixcr') to materialize an abundance PColumn with
axes [sampleId, clonotypeKey] + isAnchor annotation.

Run manually when fixtures need regeneration:

    cd blocks/titeseq-analysis/test
    uv run --with numpy python fixtures/generate_fixtures.py

Commit regenerated TSVs + manifest.json alongside the code change that
required the update.

Three variants:

- bin_mode: 5 clonotypes (Good / Partial / Failed mix) x 7 concs x 4 bins.
- no_bin_mode: 3 test clonotypes + 1 filler x 6 non-zero concs x 1 bin.
- antigen: 2 targets + 2 distractors for R4 filter verification.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

FIXTURES_DIR = Path(__file__).resolve().parent
DATA_DIR = FIXTURES_DIR / "data"
MASTER_SEED = 20260417

CONCENTRATIONS = [0.0, 1e-10, 1e-9, 3e-9, 1e-8, 1e-7, 1e-6]
BINS = [1, 2, 3, 4]

BASELINE = 1.5
AMPLITUDE = math.log(2.0)
DEFAULT_SIGMA = 0.35

TSV_HEADER = ["readCount", "nSeqCDR3", "aaSeqCDR3", "bestVGene", "bestJGene"]
V_GENE = "IGHV3-23*00"
J_GENE = "IGHJ4*00"

# Stable nt / aa alphabets used to derive per-clonotype CDR3 sequences from the
# label hash. Hashing keeps sequences unique-per-label (so MiXCR-format import
# produces distinct clonotypeKey values) without requiring manual curation.
NT_ALPHABET = "ACGT"
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def _cdr3_sequences(label: str) -> tuple[str, str]:
    # sha512 = 64 bytes, enough for 45 nt + 15 aa without cycling.
    digest = hashlib.sha512(label.encode()).digest()
    # 45 nt CDR3 = 15 AA: stays in a realistic IGHV CDR3 length range.
    nt_chars = [NT_ALPHABET[digest[i] % 4] for i in range(45)]
    nt_seq = "TGT" + "".join(nt_chars[3:])  # TGT = Cys start codon
    aa_chars = [AA_ALPHABET[digest[45 + i] % len(AA_ALPHABET)] for i in range(15)]
    aa_seq = "C" + "".join(aa_chars[1:]) + "W"
    return nt_seq, aa_seq


def _rng(label: str) -> np.random.Generator:
    seed = (MASTER_SEED + sum(ord(c) * (i + 1) for i, c in enumerate(label))) & 0xFFFFFFFF
    return np.random.default_rng(seed)


def _hill(x: float, kd: float, n: float) -> float:
    top = BASELINE + math.exp(AMPLITUDE)
    if x <= 0:
        return BASELINE
    return BASELINE + (top - BASELINE) * (x**n) / (kd**n + x**n)


def _poisson_bin_reads(
    target_bin: float,
    per_conc: int,
    sigma: float,
    rng: np.random.Generator,
) -> np.ndarray:
    centers = np.array(BINS, dtype=float)
    probs = np.exp(-0.5 * ((centers - target_bin) / sigma) ** 2)
    probs /= probs.sum()
    mean_reads = per_conc * probs
    return rng.poisson(mean_reads).astype(int)


def _write_mixcr_tsv(path: Path, rows: list[dict]) -> None:
    """Write per-sample MiXCR TSV. All rows are emitted including zero-count
    ones: the titeseq workflow's tsvFileBuilder does an outer join and emits
    empty cells for missing (sample, clonotype) pairs, which makes polars
    auto-infer the reads column as String and breaks numeric division in the
    Python pipeline. Keeping zero rows guarantees import-vdj-data's abundance
    PColumn covers every pair. Real MiXCR TSVs never contain readCount=0, but
    import-vdj-data's parser tolerates them (they simply become 0-abundance
    entries in the PColumn)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(TSV_HEADER) + "\n")
        for r in rows:
            f.write(
                f"{r['readCount']}\t{r['nSeqCDR3']}\t{r['aaSeqCDR3']}\t{V_GENE}\t{J_GENE}\n"
            )


def _build_bin_mode() -> dict:
    """Bin-mode fixture. 5 clonotypes across Good / Partial / Failed classes."""
    clonotypes = [
        {
            "label": "G_MID",
            "kd": 1e-8,
            "n": 1.0,
            "per_conc": 500,
            "sigma": DEFAULT_SIGMA,
            "expected_class": "Good",
            "expected_reason": None,
            "kd_range": [3e-9, 3e-8],
            "hill_range": [0.5, 2.0],
            "kd_out_of_range": False,
            "targets": None,
        },
        {
            "label": "G_LOW",
            "kd": 1e-9,
            "n": 1.0,
            "per_conc": 500,
            "sigma": DEFAULT_SIGMA,
            "expected_class": "Good",
            "expected_reason": None,
            "kd_range": [3e-10, 3e-9],
            "hill_range": [0.5, 2.0],
            "kd_out_of_range": False,
            "targets": None,
        },
        {
            "label": "P_NOISY",
            "kd": 1e-8,
            "n": 1.0,
            "per_conc": 60,
            "sigma": 0.85,
            "expected_class": "Partial",
            "expected_reason": None,
            "kd_range": [1e-9, 1e-7],
            "hill_range": [0.3, 10.0],
            "kd_out_of_range": False,
            "targets": None,
        },
        {
            "label": "F_HOOK",
            "kd": None,
            "n": None,
            "per_conc": 500,
            "sigma": DEFAULT_SIGMA,
            "expected_class": "Failed",
            "expected_reason": "non_monotonic_signal",
            "kd_range": None,
            "hill_range": None,
            "kd_out_of_range": None,
            # Elevated top-2 / top-3, dropped top-1 → R9b hook flag.
            "targets": [1.5, 1.5, 1.6, 1.8, 2.2, 2.8, 3.3, 3.6, 3.6, 2.5],
        },
        {
            "label": "F_INSUF_P",
            "kd": 1e-8,
            "n": 1.0,
            # Only top-3 concs pass the read floor → < min_concentration_points=5.
            "per_conc_list": [200, 1, 1, 1, 1, 1, 1, 300, 300, 300],
            "sigma": DEFAULT_SIGMA,
            "expected_class": "Failed",
            "expected_reason": "insufficient_points",
            "kd_range": None,
            "hill_range": None,
            "kd_out_of_range": None,
            "targets": None,
        },
    ]

    samples: dict[str, dict[str, int]] = {}  # sampleId → {clonotype_label: reads}
    for spec in clonotypes:
        label = spec["label"]
        rng = _rng(label)
        per_conc_list = spec.get("per_conc_list") or [spec["per_conc"]] * len(CONCENTRATIONS)
        targets = spec["targets"] or [_hill(c, spec["kd"], spec["n"]) for c in CONCENTRATIONS]
        for i, _c in enumerate(CONCENTRATIONS):
            counts = _poisson_bin_reads(targets[i], per_conc_list[i], spec["sigma"], rng)
            for j, b in enumerate(BINS):
                sample_id = f"c{i}_b{b}"
                samples.setdefault(sample_id, {})[label] = int(counts[j])

    variant_dir = DATA_DIR / "bin_mode"
    for sample_id, clone_reads in samples.items():
        rows = []
        for label, reads in clone_reads.items():
            nt, aa = _cdr3_sequences(label)
            rows.append(
                {"readCount": reads, "nSeqCDR3": nt, "aaSeqCDR3": aa, "label": label}
            )
        _write_mixcr_tsv(variant_dir / f"sample_{sample_id}.tsv", rows)

    # Per-sample metadata: concentration (Double, M) + bin (Long).
    concentration_data: dict[str, float] = {}
    bin_data: dict[str, int] = {}
    for i, c in enumerate(CONCENTRATIONS):
        for b in BINS:
            sample_id = f"c{i}_b{b}"
            concentration_data[sample_id] = float(c)
            bin_data[sample_id] = int(b)

    expected = {
        c["label"]: {
            k: c[k]
            for k in (
                "expected_class",
                "expected_reason",
                "kd_range",
                "hill_range",
                "kd_out_of_range",
            )
        }
        for c in clonotypes
    }
    # CDR3 → clonotype label lookup so tests can map PColumn rows back to
    # manifest entries without relying on clonotypeKey internals.
    cdr3_lookup = {_cdr3_sequences(c["label"])[0]: c["label"] for c in clonotypes}

    return {
        "variant": "bin_mode",
        "sample_ids": sorted(samples.keys()),
        "concentration": concentration_data,
        "bin": bin_data,
        "clonotypes": expected,
        "cdr3_to_label": cdr3_lookup,
    }


def _build_no_bin_mode() -> dict:
    """No-bin fixture. One bin per sample (so bin metadata is absent), non-zero concs only."""
    non_zero_concs = [c for c in CONCENTRATIONS if c > 0]
    n_nz = len(non_zero_concs)
    per_sample_depth = 20_000

    clonotypes = [
        {
            "label": "NB_GOOD",
            "freqs": [0.012, 0.018, 0.03, 0.05, 0.09, 0.15, 0.22, 0.25, 0.27],
            "expected_class": "Good",
            "expected_reason": None,
            "kd_range": [1e-9, 2e-7],
            "hill_range": [0.3, 2.5],
            "kd_out_of_range": False,
        },
        {
            "label": "NB_HOOK",
            "freqs": [0.012, 0.018, 0.03, 0.05, 0.09, 0.15, 0.22, 0.25, 0.08],
            "expected_class": "Failed",
            "expected_reason": "non_monotonic_signal",
            "kd_range": None,
            "hill_range": None,
            "kd_out_of_range": None,
        },
        {
            "label": "NB_LOW_R2",
            "freqs": [0.05, 0.22, 0.03, 0.20, 0.04, 0.18, 0.05, 0.17, 0.06],
            "expected_class": "Failed",
            "expected_reason_in": ["low_r2", "convergence_failure"],
            "kd_range": None,
            "hill_range": None,
            "kd_out_of_range": None,
        },
        {
            "label": "NB_FILLER",
            "freqs": None,  # absorbs remainder
        },
    ]

    rng = _rng("no_bin_mode")
    sample_ids: list[str] = []
    samples: dict[str, dict[str, int]] = {}
    for i, _c in enumerate(non_zero_concs):
        sample_id = f"nb_c{i + 1}"
        sample_ids.append(sample_id)
        row: dict[str, int] = {}
        used = 0
        for spec in clonotypes:
            if spec["label"] == "NB_FILLER":
                continue
            reads = int(rng.poisson(spec["freqs"][i] * per_sample_depth))
            row[spec["label"]] = reads
            used += reads
        row["NB_FILLER"] = max(per_sample_depth - used, 0)
        samples[sample_id] = row

    variant_dir = DATA_DIR / "no_bin_mode"
    for sample_id, clone_reads in samples.items():
        rows = []
        for label, reads in clone_reads.items():
            nt, aa = _cdr3_sequences(label)
            rows.append(
                {"readCount": reads, "nSeqCDR3": nt, "aaSeqCDR3": aa, "label": label}
            )
        _write_mixcr_tsv(variant_dir / f"sample_{sample_id}.tsv", rows)

    concentration_data = {f"nb_c{i + 1}": float(c) for i, c in enumerate(non_zero_concs)}

    expected = {
        c["label"]: {
            k: c.get(k)
            for k in (
                "expected_class",
                "expected_reason",
                "expected_reason_in",
                "kd_range",
                "hill_range",
                "kd_out_of_range",
            )
            if c.get(k) is not None
        }
        for c in clonotypes
        if c["label"] != "NB_FILLER"
    }
    cdr3_lookup = {_cdr3_sequences(c["label"])[0]: c["label"] for c in clonotypes}

    return {
        "variant": "no_bin_mode",
        "sample_ids": sample_ids,
        "concentration": concentration_data,
        "clonotypes": expected,
        "cdr3_to_label": cdr3_lookup,
    }


def _build_antigen() -> dict:
    """R4 antigen-filter fixture. Two targets + two distractors."""
    clonotypes = [
        {"label": "ANT_T1", "kd": 5e-9, "antigen": "target"},
        {"label": "ANT_T2", "kd": 5e-8, "antigen": "target"},
        {"label": "ANT_D1", "kd": 5e-9, "antigen": "other"},
        {"label": "ANT_D2", "kd": 5e-8, "antigen": "other"},
    ]

    samples: dict[str, dict[str, int]] = {}  # sampleId → {label: reads}
    sample_antigens: dict[str, str] = {}

    for spec in clonotypes:
        label = spec["label"]
        rng = _rng(label)
        targets = [_hill(c, spec["kd"], 1.0) for c in CONCENTRATIONS]
        for i, _c in enumerate(CONCENTRATIONS):
            counts = _poisson_bin_reads(targets[i], 500, 0.25, rng)
            for j, b in enumerate(BINS):
                # Samples are per (antigen, conc, bin) so antigen metadata
                # attaches cleanly to the sample axis.
                sample_id = f"ant_{spec['antigen']}_c{i}_b{b}"
                samples.setdefault(sample_id, {})[label] = int(counts[j])
                sample_antigens[sample_id] = spec["antigen"]

    variant_dir = DATA_DIR / "antigen"
    for sample_id, clone_reads in samples.items():
        rows = []
        for label, reads in clone_reads.items():
            nt, aa = _cdr3_sequences(label)
            rows.append(
                {"readCount": reads, "nSeqCDR3": nt, "aaSeqCDR3": aa, "label": label}
            )
        _write_mixcr_tsv(variant_dir / f"sample_{sample_id}.tsv", rows)

    concentration_data: dict[str, float] = {}
    bin_data: dict[str, int] = {}
    for sample_id in samples:
        # Parse the sample_id back into (conc_idx, bin) — stable format.
        _, _, c_part, b_part = sample_id.split("_")
        concentration_data[sample_id] = float(CONCENTRATIONS[int(c_part[1:])])
        bin_data[sample_id] = int(b_part[1:])

    expected = {
        "ANT_T1": {
            "expected_class": "Good",
            "kd_range": [1.5e-9, 1.5e-8],
            "hill_range": [0.5, 2.0],
            "kd_out_of_range": False,
            "antigen": "target",
        },
        "ANT_T2": {
            "expected_class": "Good",
            "kd_range": [1.5e-8, 1.5e-7],
            "hill_range": [0.5, 2.0],
            "kd_out_of_range": False,
            "antigen": "target",
        },
    }
    cdr3_lookup = {_cdr3_sequences(c["label"])[0]: c["label"] for c in clonotypes}

    return {
        "variant": "antigen",
        "sample_ids": sorted(samples.keys()),
        "concentration": concentration_data,
        "bin": bin_data,
        "antigen": sample_antigens,
        "target_antigen": "target",
        "distractor_clonotypes": ["ANT_D1", "ANT_D2"],
        "clonotypes": expected,
        "cdr3_to_label": cdr3_lookup,
    }


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    variants = {
        "bin_mode": _build_bin_mode(),
        "no_bin_mode": _build_no_bin_mode(),
        "antigen": _build_antigen(),
    }

    manifest = {
        "fixture_version": "1.0.0",
        "master_seed": MASTER_SEED,
        "concentrations_m": CONCENTRATIONS,
        "bins": BINS,
        "mixcr_tsv_headers": TSV_HEADER,
        "v_gene": V_GENE,
        "j_gene": J_GENE,
        "variants": variants,
    }

    (FIXTURES_DIR / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )

    for name, variant in variants.items():
        n_samples = len(variant["sample_ids"])
        variant_dir = DATA_DIR / name
        tsvs = sorted(p.name for p in variant_dir.glob("*.tsv"))
        print(f"  {name:<12s}: {n_samples} samples, {len(tsvs)} tsvs")


if __name__ == "__main__":
    main()
