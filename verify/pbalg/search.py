"""Search for the critical constants of an algebra.

For a gauge ``D`` the critical constants of ``A`` are

    kappa(A)  = sup_a (||a||^2 - ||a^2||) / D(a),
    lambda(A) = sup_a (|| ||a^2||1 - a^2 || - ||a^2||)_+ / D(a),

with the conventions ``0/0 = 0`` and ``t/0 = +inf`` for ``t > 0``.  ``A`` is a
PB-algebra with constant ``C`` exactly when ``max(kappa, lambda) <= C``.

Both suprema are over the (compact) unit sphere of ``A``, by homogeneity, but
the objective is nonsmooth and nonconvex.  What follows is therefore a *search*:
the number reported is a lower bound for the true supremum, obtained from many
random restarts plus a list of structured seeds, and refined until it stops
moving.  A reported value of ``0.25`` means "the search found nothing above
0.25", not "the supremum is proved to be 0.25".
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .algebras import Algebra
from .defects import pb1_defect, pb2_defect
from .gauges import (
    D1,
    Dinf,
    Dq,
    _from_coeffs,
    _orthonormal_basis,
    dist_to_scalars,
    opnorm,
)

#: Algebras for which the identity ``D_1(a) = 4 dist(a, R1)^2`` is used as a
#: fast path.  It is checked against the general search in
#: ``verify/tests/test_stampfli.py``; it is *false* for triangular algebras.
STAMPFLI_FAST_PATH = {"M2(R)", "M3(R)", "M4(R)", "H"}

UNFUNDED = float("inf")


@dataclass
class CriticalResult:
    """Outcome of a critical-constant search."""

    algebra: str
    gauge: str
    kappa: float
    lam: float
    kappa_witness: np.ndarray | None = None
    lam_witness: np.ndarray | None = None
    evaluations: int = 0

    @property
    def constant(self) -> float:
        return max(self.kappa, self.lam)

    @property
    def is_pb(self) -> bool:
        return bool(np.isfinite(self.constant))

    def summary(self) -> str:
        c = self.constant
        tag = f"C = {c:.6f}" if np.isfinite(c) else "NOT a PB-algebra"
        return (f"{self.algebra:16s} [{self.gauge}]  "
                f"kappa={self.kappa:.6f}  lambda={self.lam:.6f}   {tag}")


def make_gauge(name: str, alg: Algebra, fast: bool = True, accurate: bool = False):
    """Return a callable ``a -> D(a)`` for the named gauge on ``alg``.

    With ``accurate=False`` the ``D_1`` path skips the polish stage: cheap, but
    it can *under*estimate the budget on a proper subalgebra and so
    *over*estimate a ratio.  Every coarse search is therefore followed by a
    refinement pass with ``accurate=True`` (see :func:`critical_constants`).
    """
    if name == "D1":
        if fast and alg.name in STAMPFLI_FAST_PATH:
            return lambda a, rng=None: 4.0 * dist_to_scalars(a) ** 2
        if accurate:
            # The polish stage dominates the cost and grows with dim(A): each
            # restart is a derivative-free run over dim(A) parameters.  Near the
            # optimum a couple of restarts suffice, so spend them where they are
            # cheap.
            polish = 6 if alg.dim <= 4 else (3 if alg.dim <= 6 else 2)
            return lambda a, rng=None: D1(alg, a, restarts=10, iters=80,
                                          polish=polish, rng=rng)
        return lambda a, rng=None: D1(alg, a, restarts=6, iters=60, polish=0,
                                      rng=rng)
    if name == "Dinf":
        return lambda a, rng=None: Dinf(alg, a, restarts=8, rng=rng)
    if name.startswith("D") and name[1:].isdigit():
        q = int(name[1:])
        return lambda a, rng=None: Dq(alg, a, q, restarts=8, rng=rng)
    raise ValueError(f"unknown gauge {name!r}")


def _seeds(alg: Algebra, rng: np.random.Generator, tries: int) -> list[np.ndarray]:
    """Random coefficient vectors plus the standard structured witnesses."""
    q = _orthonormal_basis(alg)
    dim = q.shape[1]
    seeds = [rng.standard_normal(dim) for _ in range(tries)]
    for b in alg.basis:  # each basis element
        seeds.append(q.T @ b.ravel())
    n = alg.size
    for i in range(n - 1):  # matrix units and Jordan blocks
        e = np.zeros((n, n))
        e[i, i + 1] = 1.0
        if alg.contains(e):
            seeds.append(q.T @ e.ravel())
            seeds.append(q.T @ (np.eye(n) + e).ravel())
    return seeds


def _sup_ratio(alg: Algebra, defect, gauge, rng: np.random.Generator, tries: int,
               maxiter: int) -> tuple[float, np.ndarray | None, int]:
    q = _orthonormal_basis(alg)
    n = alg.size
    best: float = 0.0
    witness: np.ndarray | None = None
    calls = 0

    def ratio(c: np.ndarray) -> float:
        nonlocal calls
        calls += 1
        a = _from_coeffs(q, c, n)
        na = opnorm(a)
        if na < 1e-10:
            return 0.0
        a = a / na
        d = defect(a)
        if d <= 1e-10:
            return 0.0
        g = gauge(a, rng)
        if g <= 1e-9:
            return 1e9  # a defect with no budget: not a PB-algebra at any C
        return d / g

    for c0 in _seeds(alg, rng, tries):
        res = minimize(lambda c: -ratio(c), c0, method="Nelder-Mead",
                       options={"maxiter": maxiter, "xatol": 1e-9, "fatol": 1e-11})
        val = -float(res.fun)
        if val > best:
            a = _from_coeffs(q, res.x, n)
            na = opnorm(a)
            best, witness = val, (a / na if na > 1e-10 else None)
        if best >= 1e9:
            return UNFUNDED, witness, calls
    return best, witness, calls


def _refine(alg: Algebra, defect, witness: np.ndarray | None, gauge_acc,
            rng: np.random.Generator, maxiter: int) -> tuple[float, np.ndarray | None]:
    """Re-evaluate and locally re-maximise a ratio using the accurate gauge."""
    if witness is None:
        return 0.0, None
    q = _orthonormal_basis(alg)
    n = alg.size

    def ratio(c: np.ndarray) -> float:
        a = _from_coeffs(q, c, n)
        na = opnorm(a)
        if na < 1e-10:
            return 0.0
        a = a / na
        d = defect(a)
        if d <= 1e-10:
            return 0.0
        g = gauge_acc(a, rng)
        if g <= 1e-9:
            return 1e9
        return d / g

    c0 = q.T @ witness.ravel()
    res = minimize(lambda c: -ratio(c), c0, method="Nelder-Mead",
                   options={"maxiter": maxiter, "xatol": 1e-9, "fatol": 1e-11})
    best = max(ratio(c0), -float(res.fun))
    a = _from_coeffs(q, res.x, n)
    na = opnorm(a)
    return best, (a / na if na > 1e-10 else witness)


def critical_constants(alg: Algebra, gauge: str = "D1", tries: int = 30,
                       maxiter: int = 400, seed: int = 0, fast: bool = True,
                       refine_iter: int = 120) -> CriticalResult:
    """Search for ``kappa`` and ``lambda`` of ``alg`` against ``gauge``.

    Two passes: a broad coarse search, then a local refinement of the winning
    witness with the accurate gauge.  The refinement matters -- a coarse gauge
    that undershoots the budget inflates the ratio, which is how an earlier run
    of this code produced a spurious ``kappa(T2(R)) = 0.317``.
    """
    rng = np.random.default_rng(seed)
    g = make_gauge(gauge, alg, fast=fast)
    g_acc = make_gauge(gauge, alg, fast=fast, accurate=True)

    kap, kw, c1 = _sup_ratio(alg, pb1_defect, g, rng, tries, maxiter)
    lam, lw, c2 = _sup_ratio(alg, pb2_defect, g, rng, tries, maxiter)

    # Each refinement step costs an accurate gauge evaluation, which is itself a
    # maximisation over dim(A) parameters; without this the refinement of a
    # six-dimensional algebra takes longer than the rest of the run together.
    if alg.dim > 4 and alg.name not in STAMPFLI_FAST_PATH:
        refine_iter = min(refine_iter, 60 if alg.dim <= 6 else 40)

    if np.isfinite(kap):
        kap, kw = _refine(alg, pb1_defect, kw, g_acc, rng, refine_iter)
    if np.isfinite(lam):
        lam, lw = _refine(alg, pb2_defect, lw, g_acc, rng, refine_iter)

    return CriticalResult(alg.name, gauge, kap, lam, kw, lw, c1 + c2)
