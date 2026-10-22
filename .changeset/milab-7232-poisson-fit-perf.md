---
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

Faster bin-mode fits with less memory; the results are unchanged

- The Adams 2016 fit runs 1.6–2x faster and needs 4–6x less peak memory per worker process.
- Fit worker processes follow the job's CPU limit instead of the node's core count, for both curve models.
