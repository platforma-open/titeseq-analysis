# CID Conflict Investigation — titeseq-analysis

**Date:** 2026-04-24
**Symptom:** Occasional "CID conflict" errors when running the titeseq-analysis block.

---

## 1. What CID actually means

CID = **Canonical ID** — Platforma's content-addressable hash for a resource. The backend throws `CIDConflictError` (`core/pl/platform/core/transaction/resource_error.go:403`) when a field is populated with one CID in one transaction attempt and a different CID in a later attempt:

> `"CID conflict in field %q (resource type %q). Current field CID: %q, new field CID: %q"`

The backend has a self-recovery mechanism (`field_state_change.go:307-313`, "zebra CID conflicts"), which is why the errors are **transient/occasional** rather than persistent.

**CID conflicts are caused by non-deterministic resource content OR by the same logical field receiving two different resources across parallel/retry paths.**

---

## 2. Root cause in titeseq-analysis

### The dual-import pattern

`workflow/src/main.tpl.tengo` imports each of `mean_bin.tsv` and `fitted_mean_bin.tsv` **twice** — once with 3 axes for the internal view, once with 2 axes for the export view:

```tengo
meanBinInternal := importSignal(fit.getFile("mean_bin.tsv"), "meanBin", meanBinInternalColSpec, internalConcAxes)   // 3 axes
meanBinExport   := importSignal(fit.getFile("mean_bin.tsv"), "meanBin", meanBinExportColSpec,   exportConcAxes)     // 2 axes
fittedMeanBinInternal := importSignal(fit.getFile("fitted_mean_bin.tsv"), "fittedMeanBin", fittedMeanBinInternalColSpec, internalConcAxes)
fittedMeanBinExport   := importSignal(fit.getFile("fitted_mean_bin.tsv"), "fittedMeanBin", fittedMeanBinExportColSpec,   exportConcAxes)
```

Both variants share:
- `name`: `pl7.app/vdj/meanBin` (or `pl7.app/vdj/fittedMeanBin`)
- `valueType`: `Double`
- `domain`: `{"pl7.app/block": blockId}`

They differ only in:
- `axesSpec` length (3 vs 2)
- One annotation (`pl7.app/hideDataFromGraphs` on the export variant)

### The comment in main.tpl.tengo already flags the risk

> "Export variants feed exports.pf (2-axis view for downstream discovery) and carry `pl7.app/hideDataFromGraphs` so `createPFrameForGraphs` excludes them from this block's own pframe — otherwise **both copies share name/valueType/domain/annotations**, and findColumnBy resolves the y-source to the 2-axis export sibling, silently dropping the numeric-concentration x default."

The `hideDataFromGraphs` annotation fixes UI-level pcolumn resolution, but it does **not** change the underlying resource identity. The two variants remain two distinct resources advertising the same `(name, valueType, domain)` tuple to the result pool.

### Why the conflict is intermittent

Both variants render through `xsv.importFile(...)` on the same input file. When the backend's deduplication/merge layer sees two pcolumn resources with identical logical identity but different resource CIDs being added to **exported** pframes (both `signalPf` and `exportPf` are surfaced via `pframes.exportFrame`/`exports.pf`), a downstream consumer's field can be populated by whichever path resolves first — producing CID drift across retries. The backend's self-recovery kicks in, so most attempts eventually succeed — hence "occasional."

### Why Python output determinism is NOT the bug

`software/src/pipeline.py` sorts outputs deterministically (`sorted_points = fit_points.sort([COL_CLONOTYPE, COL_CONC_VAL])`, `clonotypes_sorted = sorted(all_clonotypes)`). The TSV content is stable.

---

## 3. Idiomatic patterns in other blocks

### `blocks/mixcr-clonotyping/workflow/src/calculate-export-specs.lib.tengo`

Deterministic domain canonicalization:

```tengo
toCombinedDomainValue := func(spec) {
    result := [spec.name]
    // getKeys sort keys
    for domain in maps.getKeys(spec.domain) {
        result = append(result, [domain, spec.domain[domain]])
    }
    return result
}
```

`maps.getKeys` internally `slices.quickSortInPlace`s — never iterate a map directly if the iteration feeds canonical content.

### `blocks/clonotype-clustering/workflow/src/main.tpl.tengo`

**One import per output file.** Each logical output (`clusterToSeqPf`, `clusterRadiusPf`, `abundancesPerClusterPf`) is imported exactly once from its own TSV. When columns need to appear in multiple pFrames (`opf`, `bubblePlotPf`, `msaPf`, `epf`), the SAME pcolumn resource is reused — no duplicate imports, no duplicate identities.

### Standard guidance

- Each logical output value = one TSV file from the Python/tool layer
- Each TSV = exactly one `xsv.importFile` call
- Multiple views/pframes reuse the resulting pcolumn resource; never re-import with different axes

---

## 4. Recommended fix (ranked)

### A. Differentiate the column identity (lowest risk, smallest diff)

Give the two variants distinct `(name, domain)` tuples so the result pool never has to choose between them.

In `main.tpl.tengo`, add a view discriminator to the export variants' domain:

```tengo
meanBinExportColSpec := {
    name: "pl7.app/vdj/meanBin",
    valueType: "Double",
    domain: {
        "pl7.app/block": blockId,
        "pl7.app/titeseq/view": "export"    // <-- new
    },
    annotations: { ...existing... }   // hideDataFromGraphs becomes redundant but harmless
}
```

…and add the matching `"pl7.app/titeseq/view": "internal"` to `meanBinInternalColSpec` (and both fitted variants).

Downstream blocks that currently pick up the export variant by `(name, domain subset)` keep working because domain matching is a subset match. Own-block Graph Maker stops seeing the export sibling as a valid y-source because the domain now disagrees.

### C. Import once, re-spec for export (most idiomatic)

Import each TSV with the 3-axis internal config **once**. Build the export view by dropping the `concentrationAM` axis via an SDK helper (if one exists — check `pSpec.setAxisSpec` / `pframes.aggregate` / axis-collapse helpers), or by emitting the export-axis view as a separate TSV from Python and importing it as a distinct logical column.

This removes the dual-import entirely and aligns with clonotype-clustering's pattern.

### D. Don't expose the internal view to the pool

If `signalPf`'s columns are only consumed by this block's own Graph Maker, see whether they can be delivered via a non-`exportFrame` output path (e.g., handle-based resolution inside `ctx.outputs.resolve`) so they never enter the cross-block result pool. This requires confirming with the platform team that such a path exists for pframes; `pframes.exportFrame` is how we've always done it.

**Recommendation:** do **A** first (single-line import swap, matches SDK convention used by 5+ other blocks). If conflicts persist, layer on **B**. **C** is the cleanest long-term refactor.

---

## 6. Preventive guardrails for future blocks

1. **Use `canonical.encode`, not `json.encode`, for anything that feeds a resource's bytes.** `writeFile`, domain values that embed config JSON, cache keys, trace payloads — all of these should be `canonical.encode(...)` when the input is a map. `json.encode` is fine for strings, numbers, arrays, and for printf-style debug output.
2. **Never re-import the same TSV with different axis configs.** One TSV → one `xsv.importFile`. If two views are needed, differentiate them at the spec/domain level or emit two TSVs from the tool.
3. **Make (name, domain) globally unique in the pool.** Use a `pl7.app/<block>/view` domain key (or similar) when you need multiple parallel presentations of the same logical quantity.
4. **Sort before you iterate.** When building content from a Tengo map, use `maps.getKeys(m)` (returns sorted keys) rather than `for k, v in m`, whenever the order affects output bytes.
5. **Review code that comments "both copies share name/valueType/domain".** That phrase in a workflow is a CID-conflict smell.

---

## Files referenced

**titeseq-analysis:**
- `blocks/titeseq-analysis/workflow/src/main.tpl.tengo:99` (`json.encode(paramsPayload)` — switch to `canonical.encode`)
- `blocks/titeseq-analysis/workflow/src/main.tpl.tengo` (dual-import structure — spec defs ~115-180, imports ~260-300, pframe builders ~320-365)

**SDK helpers:**
- `core/platforma/sdk/workflow-tengo/src/canonical.lib.tengo` — `canonical.encode` (sorted-key JSON)
- `core/platforma/sdk/workflow-tengo/src/pframes/spec.lib.tengo` — `makeTrace`, `createSpecDistiller`, `prepareAxisFiltersAndAxesSpec`
- `core/platforma/sdk/workflow-tengo/src/anonymize/index.lib.tengo` — `anonymizeFields`, `anonymizePKeys`
- `core/platforma/sdk/workflow-tengo/src/pframes/slice-data.tpl.tengo` — axis-filter view derivation
- `core/platforma/sdk/workflow-tengo/src/maps.lib.tengo` — `getKeys` (sorted)

**Reference patterns (other blocks using `canonical.encode`):**
- `blocks/repertoire-distance/workflow/src/run-distance.tpl.tengo:22` — `writeFile("metrics.json", canonical.encode(metrics))`
- `blocks/clonotype-browser/workflow/src/annotations.lib.tengo`, `export-table.tpl.tengo`
- `blocks/clonotype-enrichment/workflow/src/pf-enrichment-conv-export.lib.tengo`, `pf-frequency-conv.lib.tengo`
- `blocks/repertoire-diversity/workflow/src/main.tpl.tengo`
- `blocks/cell-browser/workflow/src/annotations.lib.tengo`

**Comparison blocks:**
- `blocks/mixcr-clonotyping/workflow/src/calculate-export-specs.lib.tengo` (deterministic domain canonicalization via `maps.getKeys`)
- `blocks/clonotype-clustering/workflow/src/main.tpl.tengo` (single-import idiom)

**Backend:**
- `core/pl/platform/core/transaction/resource_error.go:400` — `CIDConflictError`
- `core/pl/platform/core/transaction/field_state_change.go:307` — self-recovery ("zebra") mechanism
