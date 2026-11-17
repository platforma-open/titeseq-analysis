---
"@platforma-open/platforma-open.titeseq-analysis": patch
"@platforma-open/platforma-open.titeseq-analysis.model": patch
"@platforma-open/platforma-open.titeseq-analysis.ui": patch
---

Replicate panel: choose the columns, and hide the empty ones.

- The panel that opens from a variant's row has a Columns list. Use it to show or hide a column.
  The table's own Columns button is in the page title, behind the panel.
- By default the panel shows the same columns as before, less each column that has no value for
  that variant. For example, Fit failure reason is not shown when every replicate of the variant
  was fitted.
- A column that you show or hide in the list stays that way for the other variants.
