---
'@platforma-open/platforma-open.titeseq-analysis.model': patch
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

Hide metadata columns that are empty for every sample from the Antigen concentration and FACS bin dropdowns — selecting one previously bombed the Python pipeline with a confusing dtype error. When such a column is still selected (e.g. via a saved project), the bin validator now raises a friendly "pick a populated column" message with the same wording as the concentration guard.
