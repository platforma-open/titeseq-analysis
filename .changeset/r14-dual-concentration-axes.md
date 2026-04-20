---
'@platforma-open/platforma-open.titeseq-analysis.workflow': patch
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

R14 compliance: declare both canonical String axis `concentrationStr` and numeric Long axis `concentrationAM` (attomolar) on `meanBin` and `fittedMeanBin` PColumns. The canonical string is the join key (preserving upstream metadata byte-for-byte, preventing float-serialization drift); the Long sibling continues to drive Graph Maker's log-scale X rendering.
