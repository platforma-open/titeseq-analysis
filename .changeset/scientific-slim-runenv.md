---
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

Switch Python runtime to `runenv-python-3:3.12.10-scientific-slim` (bundles polars-lts-cpu, numpy, scipy, pyarrow only). Drops the shipped runtime from ~2.5 GB to ~100 MB by dropping the base runenv's TensorFlow/torch/transformers/cudf stack — none of which this block uses. The `software/src/requirements.txt` pins match the wheels bundled by the slim variant, so runtime `pip install --no-index --find-links={runenv}/packages/` resolves cleanly.
