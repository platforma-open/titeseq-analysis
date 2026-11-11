---
"@platforma-open/platforma-open.titeseq-analysis": major
"@platforma-open/platforma-open.titeseq-analysis.model": major
"@platforma-open/platforma-open.titeseq-analysis.ui": major
"@platforma-open/platforma-open.titeseq-analysis.workflow": major
"@platforma-open/platforma-open.titeseq-analysis.kind": major
"@platforma-open/platforma-open.titeseq-analysis.software": major
"@platforma-open/platforma-open.titeseq-analysis.test": major
---

Publish the fit software again under a new version. The 2.3.0 archive in the software registry
holds an older build, because that name was taken before the release and the upload was skipped.
The 2.6.0 workflow passes options that build does not know, so a run on a backend that uses the
archive stops with "unrecognized arguments".
