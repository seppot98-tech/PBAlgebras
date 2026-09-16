"""The two axiom defects and their ratios against a gauge."""

from __future__ import annotations

import numpy as np

from .gauges import opnorm


def pb1_defect(a: np.ndarray) -> float:
    """``||a||^2 - ||a^2||``, the failure of the square property."""
    return opnorm(a) ** 2 - opnorm(a @ a)


def pb2_defect(a: np.ndarray) -> float:
    """``|| ||a^2|| 1 - a^2 || - ||a^2||``, the failure of spectral reality.

    Nonpositive exactly when the unrelaxed condition (SR) holds at ``a``.
    """
    n = a.shape[0]
    sq = a @ a
    s = opnorm(sq)
    return opnorm(s * np.eye(n) - sq) - s


def pb1_ratio(a: np.ndarray, gauge_value: float, tol: float = 1e-9) -> float:
    """``(||a||^2 - ||a^2||) / D(a)``; ``inf`` if the budget is zero but the defect is not."""
    d = pb1_defect(a)
    if gauge_value <= tol:
        return 0.0 if d <= tol else float("inf")
    return d / gauge_value


def pb2_ratio(a: np.ndarray, gauge_value: float, tol: float = 1e-9) -> float:
    """``(|| ||a^2||1 - a^2 || - ||a^2||) / D(a)``, clipped below at zero."""
    d = pb2_defect(a)
    if d <= tol:
        return 0.0
    if gauge_value <= tol:
        return float("inf")
    return d / gauge_value
