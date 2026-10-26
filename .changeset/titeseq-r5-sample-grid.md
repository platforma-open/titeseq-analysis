---
'@platforma-open/platforma-open.titeseq-analysis.software': patch
'@platforma-open/platforma-open.titeseq-analysis.workflow': patch
'@platforma-open/platforma-open.titeseq-analysis': patch
---

Stop the Fit Log reporting missing bin and concentration data that is present

The "missing (bin, concentration) combinations" warning fired on every real dataset. It
counted the rows of each clonotype, and upstream output is sparse: a clonotype with no reads
in a bin has no row there. That is a true zero, it contributes nothing to the mean bin, and it
changes no fitted value.

The warning now reports what it was meant to report — a (bin, concentration) sample that was
never measured, which does drop a bin from the mean bin at that concentration. It reads the
sample list for this, because the reads table cannot distinguish an unmeasured sample from a
bin that is legitimately empty, such as the high-fluorescence bins at 0 M.
