"""Goodness of fit shared by the least-squares equations."""

from __future__ import annotations

import numpy as np


def weighted_r2(y: np.ndarray, y_hat: np.ndarray, w: np.ndarray) -> float:
    """R²_w = 1 − Σ(w·(y−ŷ)²) / Σ(w·(y−ȳ_w)²), ȳ_w = Σ(w·y)/Σ(w).

    Zero-weight points contribute to neither sum. Negative values are not clamped. Returns nan when
    Σw = 0 or when the data have no spread. Unit weights give the ordinary R².
    """
    w = np.asarray(w, dtype=float)
    y = np.asarray(y, dtype=float)
    y_hat = np.asarray(y_hat, dtype=float)
    sum_w = float(w.sum())
    if sum_w == 0.0:
        return float("nan")
    y_bar_w = float((w * y).sum() / sum_w)
    num = float((w * (y - y_hat) ** 2).sum())
    den = float((w * (y - y_bar_w) ** 2).sum())
    if den == 0.0:
        return float("nan")
    return 1.0 - num / den
