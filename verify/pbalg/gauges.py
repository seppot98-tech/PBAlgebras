"""The commutator gauges.

For an element ``a`` of an algebra ``A`` the gauges are

    D_1(a)    = sup_{||z||=1} ||[a,z]||^2            = ||ad_a||^2
    D_q(a)    = sup_{||z||=1} || [a,z]^q ||^(2/q)
    D_inf(a)  = sup_{||z||=1} r([a,z])^2

Each is homogeneous of degree two in ``a`` and invariant under translation by
scalars, and ``D_inf <= D_q <= D_1 <= 4||a||^2``.

All three are suprema of a *convex* (``D_1``) or merely continuous
(``D_q``, ``D_inf``) function over the unit ball of ``A``, so the optimisation
is nonconvex maximisation.  We use scale-invariant ratios together with many
random restarts; for ``D_1`` the restarts are refined by a Frank-Wolfe style
linearisation step, which converges quickly and is exact for the full matrix
algebra.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .algebras import Algebra


def opnorm(m: np.ndarray) -> float:
    """Operator norm (largest singular value)."""
    return float(np.linalg.norm(m, 2))


def spectral_radius(m: np.ndarray) -> float:
    return float(np.max(np.abs(np.linalg.eigvals(m))))


def _orthonormal_basis(alg: Algebra) -> np.ndarray:
    """Frobenius-orthonormal basis of the span, as columns of a matrix."""
    flat = np.stack([b.ravel() for b in alg.basis], axis=1)
    q, _ = np.linalg.qr(flat)
    return q


def _from_coeffs(q: np.ndarray, c: np.ndarray, n: int) -> np.ndarray:
    return (q @ c).reshape(n, n)


def _polar(g: np.ndarray) -> np.ndarray:
    """Orthogonal polar factor ``UV^T`` of ``g``; the linear oracle for the
    operator-norm unit ball."""
    u, _, vt = np.linalg.svd(g)
    return u @ vt


def D1(alg: Algebra, a: np.ndarray, restarts: int = 16, iters: int = 80,
       polish: int = 8, rng: np.random.Generator | None = None) -> float:
    """``D_1(a) = ||ad_a||^2``, the supremum taken over the unit ball of ``alg``.

    Maximisation of the convex function ``z -> ||[a,z]||`` over the unit ball is
    done in two stages.

    *Frank-Wolfe.*  Linearise at the current point -- the subgradient is
    ``ad_a^*(uv^T)`` for ``(u,v)`` the top singular pair of ``[a,z]`` -- and move
    to the best available extreme point.  For the full matrix algebra the linear
    oracle is exact: the maximiser of ``<g,z>`` over the operator-norm ball is
    the orthogonal polar factor of ``g``.

    *Polish.*  For a **proper subalgebra** the polar factor need not lie in the
    subalgebra, so the oracle is only approximate and the iteration can stall
    below the true supremum (it does, on ``T_n(R)``, by several percent).  We
    therefore follow it with a derivative-free maximisation of the
    scale-invariant ratio ``||[a,z]||^2 / ||z||^2``, seeded by the Frank-Wolfe
    optimum and by random restarts.  The result is the larger of the two stages,
    hence always a lower bound for ``D_1(a)``, and it agrees with brute force to
    ``1e-9`` on the test algebras.
    """
    rng = rng or np.random.default_rng(0)
    q = _orthonormal_basis(alg)
    n = alg.size
    best = 0.0
    best_z = None

    def obj(z: np.ndarray) -> float:
        return opnorm(a @ z - z @ a)

    for _ in range(restarts):
        z = _from_coeffs(q, rng.standard_normal(q.shape[1]), n)
        nz = opnorm(z)
        if nz < 1e-14:
            continue
        z = z / nz
        val = obj(z)
        for _ in range(iters):
            w = a @ z - z @ a
            u, _, vt = np.linalg.svd(w)
            g = np.outer(u[:, 0], vt[0])
            grad = a.T @ g - g @ a.T
            projected = _from_coeffs(q, q.T @ grad.ravel(), n)
            candidates = []
            for cand in (projected,
                         _from_coeffs(q, q.T @ _polar(projected).ravel(), n)):
                nc = opnorm(cand)
                if nc > 1e-14:
                    candidates.append(cand / nc)
            if not candidates:
                break
            new_z = max(candidates, key=obj)
            new_val = obj(new_z)
            if new_val <= val + 1e-12:
                break
            z, val = new_z, new_val
        if val > best:
            best, best_z = val, z

    if polish:
        def neg(c: np.ndarray) -> float:
            z = _from_coeffs(q, c, n)
            nz = opnorm(z)
            if nz < 1e-12:
                return 0.0
            return -obj(z / nz)

        seeds = []
        if best_z is not None:
            seeds.append(q.T @ best_z.ravel())
        seeds.extend(rng.standard_normal(q.shape[1]) for _ in range(polish))
        for c0 in seeds:
            res = minimize(neg, c0, method="Nelder-Mead",
                           options={"maxiter": 2000, "xatol": 1e-11,
                                    "fatol": 1e-13})
            best = max(best, -float(res.fun))

    return best ** 2


def _ratio_search(alg: Algebra, f, restarts: int, rng: np.random.Generator) -> float:
    """Maximise a scale-invariant ratio ``f(z)`` over ``z`` in ``alg``."""
    q = _orthonormal_basis(alg)
    n = alg.size
    best = 0.0
    for _ in range(restarts):
        c0 = rng.standard_normal(q.shape[1])

        def neg(c: np.ndarray) -> float:
            z = _from_coeffs(q, c, n)
            nz = opnorm(z)
            if nz < 1e-12:
                return 0.0
            return -f(z / nz)

        res = minimize(neg, c0, method="Nelder-Mead",
                       options={"maxiter": 3000, "xatol": 1e-10, "fatol": 1e-12})
        best = max(best, -float(res.fun))
    return best


def Dq(alg: Algebra, a: np.ndarray, qexp: float, restarts: int = 24,
       rng: np.random.Generator | None = None) -> float:
    """``D_q(a)`` for finite ``q >= 1``."""
    rng = rng or np.random.default_rng(0)
    k = int(round(qexp))
    if abs(k - qexp) > 1e-12 or k < 1:
        raise ValueError("Dq is implemented for integer q >= 1")

    def f(z: np.ndarray) -> float:
        w = a @ z - z @ a
        return opnorm(np.linalg.matrix_power(w, k)) ** (2.0 / k)

    return _ratio_search(alg, f, restarts, rng)


def Dinf(alg: Algebra, a: np.ndarray, restarts: int = 24,
         rng: np.random.Generator | None = None) -> float:
    """``D_inf(a) = sup_z r([a,z])^2``."""
    rng = rng or np.random.default_rng(0)

    def f(z: np.ndarray) -> float:
        return spectral_radius(a @ z - z @ a) ** 2

    return _ratio_search(alg, f, restarts, rng)


def dist_to_scalars(a: np.ndarray) -> float:
    """``dist(a, R1)`` in the operator norm."""
    n = a.shape[0]
    r = opnorm(a)
    res = minimize_scalar(
        lambda t: opnorm(a - t * np.eye(n)),
        bounds=(-r - 1.0, r + 1.0),
        method="bounded",
        options={"xatol": 1e-13},
    )
    return float(res.fun)
