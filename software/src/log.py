"""Timestamped stdout progress logging for the Fit Log UI.

Writes to stdout with `flush=True` so the Tengo `saveStdoutStream()` output
streams live into `PlLogView` (consumed by `app.model.outputs.logHandle`).
Elapsed seconds are measured from the first `log()` call in the process — the
pipeline is typically bracketed by a single CLI invocation, so elapsed time
reads as "seconds since Python started doing work".
"""

from __future__ import annotations

import sys
import time

_START: float | None = None


def log(message: str) -> None:
    global _START
    now = time.monotonic()
    if _START is None:
        _START = now
    elapsed = now - _START
    print(f"[{elapsed:7.2f}s] {message}", flush=True, file=sys.stdout)
