---
"@platforma-open/platforma-open.titeseq-analysis.software": patch
"@platforma-open/platforma-open.titeseq-analysis.workflow": patch
"@platforma-open/platforma-open.titeseq-analysis": patch
---

Publish the fit software again under a new version. The 2.3.0 archive in the software registry
holds an older build, because that name was taken before the release and the upload was skipped.
The 2.6.0 workflow passes options that build does not know, so a run on a backend that uses the
archive stops with "unrecognized arguments".
