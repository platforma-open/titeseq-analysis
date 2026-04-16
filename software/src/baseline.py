"""R6: global baseline B = arithmetic mean of per-clonotype mean_bin at c=0.

Clonotypes whose c=0 point is filtered by the read floor are excluded from the average.
Returns None when no c=0 data survives (downstream uses 4-param fit).
"""

from __future__ import annotations

import polars as pl

from normalization import SIGNAL


def compute_global_baseline(c0_points: pl.DataFrame) -> float | None:
    """Arithmetic mean of signal at c=0 across clonotypes that survived R8 floor.

    c0_points is the floor-passed subset of the signal frame where concentration == 0.
    Returns None if the frame is empty.
    """
    if c0_points.height == 0:
        return None
    mean_val = c0_points.select(pl.col(SIGNAL).mean()).item()
    if mean_val is None:
        return None
    return float(mean_val)
