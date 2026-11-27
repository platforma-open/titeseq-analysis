---
'@platforma-open/platforma-open.titeseq-analysis': minor
'@platforma-open/platforma-open.titeseq-analysis.model': minor
'@platforma-open/platforma-open.titeseq-analysis.ui': minor
'@platforma-open/platforma-open.titeseq-analysis.workflow': minor
'@platforma-open/platforma-open.titeseq-analysis.software': minor
---

- New "Run Quality" page: fit outcomes, failed fits and raw counts, per replicate and for all replicates together
- The variant panel opens for every run, with a "Curves" tab that shows the variant's fitted titration curves, Kd,app and its standard error
- The "Replicates" table hides the columns that are empty for the variant and adds log10 Kd,app per replicate
- The tables show log10 Kd,app by default; -log10 Kd,app and "Kd,app out of range" are optional columns
- A variant whose replicates or synonymous variants all failed for one reason names that reason
- Removed the "Kd Distribution" and "Titration Curves" pages
- Small bug fixes and improvements
