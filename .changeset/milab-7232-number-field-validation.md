---
"@platforma-open/platforma-open.titeseq-analysis": patch
"@platforma-open/platforma-open.titeseq-analysis.model": patch
"@platforma-open/platforma-open.titeseq-analysis.ui": patch
---

Every number field in the settings drawer states each rule where the field itself checks it: its
bounds are bounds, and a custom check covers only what a bound cannot say — a whole number, and the
two paired limits (R² Good against Failed, Hill coefficient min against max). The bound messages
come from the field component, so they read as they do in every other block.

Run is disabled while the plateau floor is empty, as it already was for every other fit setting.

A concentration column label with spaces no longer raises a warning. The label still becomes the
Kd,app unit; spaces in it break nothing.

Clearing a number field leaves it empty and marked "Value is required", and disables Run until a
value is entered. A cleared box used to read empty over the number the run still used.
