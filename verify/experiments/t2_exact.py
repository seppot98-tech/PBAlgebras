"""An exact critical-constant computation for T_2(R).

Whether ``T_2(R)`` is a PB-algebra at ``C = 1/4`` is a headline claim -- it
decides whether any constant separates ``T_2(R)`` from ``M_2(R)`` -- so it
deserves better than a generic nonconvex search.

Every commutator in ``T_2(R)`` is a multiple of ``E_12``: for
``x = [[al, ga], [0, be]]`` and ``z = [[p, c], [0, q]]``,

    [x, z] = ( c(al - be) + ga(q - p) ) E_12 ,

so ``||[x,z]||`` is the absolute value of a *linear* functional of ``z``, and

    D_1(x) = S(al - be, ga)^2,   S(s,t) = max{ |s c + t(q-p)| : ||z|| <= 1 } .

``S`` is the support function of a convex body, so computing it is maximisation
of a linear functional over a convex set: no local maxima.  Better, it has a
closed dual form.  Since ``z_21 = 0``, the functional is ``<G, z>`` for *any*
``G`` whose (1,1), (1,2), (2,2) entries are ``-t, s, t`` -- the (2,1) entry is
free -- and duality for the operator-norm ball gives

    S(s, t) = min_{g in R} || [[-t, s], [g, t]] ||_*        (nuclear norm),

a one-dimensional convex minimisation.  This is exact to machine precision and
costs microseconds, which makes a dense sweep over ``x`` feasible.

Run directly to reproduce the numbers quoted in the paper.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pbalg.gauges import opnorm  # noqa: E402


def support(s: float, t: float) -> float:
    """``max{ |s c + t(q - p)| : || [[p,c],[0,q]] || <= 1 }``, by duality."""
    def nuclear(g: float) -> float:
        m = np.array([[-t, s], [g, t]])
        return float(np.linalg.svd(m, compute_uv=False).sum())

    res = minimize_scalar(nuclear, bounds=(-50.0, 50.0), method="bounded",
                          options={"xatol": 1e-13})
    return float(res.fun)


def D1_exact(x: np.ndarray) -> float:
    """``D_1(x)`` for ``x`` in ``T_2(R)``, exactly."""
    return support(x[0, 0] - x[1, 1], x[0, 1]) ** 2


def D1_primal(x: np.ndarray, samples: int = 20000,
              seed: int = 0) -> float:
    """Independent check on :func:`D1_exact` by sampling the primal problem."""
    rng = np.random.default_rng(seed)
    s, t = x[0, 0] - x[1, 1], x[0, 1]
    best = 0.0
    for _ in range(samples):
        p, c, q = rng.standard_normal(3)
        z = np.array([[p, c], [0.0, q]])
        n = opnorm(z)
        if n < 1e-12:
            continue
        p, c, q = p / n, c / n, q / n
        best = max(best, abs(s * c + t * (q - p)))
    return best ** 2


def defects(x: np.ndarray) -> tuple[float, float]:
    n2 = opnorm(x) ** 2
    sq = x @ x
    s = opnorm(sq)
    return n2 - s, opnorm(s * np.eye(2) - sq) - s


def sweep(grid: int = 601, span: float = 4.0) -> dict:
    """Maximise both ratios over the unit sphere of ``T_2(R)``.

    After normalisation the problem has two parameters, not three.  Elements
    with ``ga = 0`` are diagonal, and there both defects are nonpositive, so
    they contribute nothing; every other element is a scalar multiple of one
    with ``ga = 1``, and the ratios are scale invariant.  So it suffices to
    sweep ``(al, be)`` with ``ga = 1``.
    """
    best: dict[str, tuple[float, np.ndarray | None]] = {
        "kappa": (0.0, None), "lambda": (0.0, None)
    }
    axis = np.linspace(-span, span, grid)
    for al in axis:
        for be in axis:
            x = np.array([[al, 1.0], [0.0, be]])
            x = x / opnorm(x)
            d1, d2 = defects(x)
            if d1 <= 1e-12 and d2 <= 1e-12:
                continue
            g = D1_exact(x)
            for key, d in (("kappa", d1), ("lambda", d2)):
                if d <= 1e-12:
                    continue
                r = d / g if g > 1e-12 else float("inf")
                if r > best[key][0]:
                    best[key] = (r, x.copy())
    return best


def refine(x0: np.ndarray, which: str) -> tuple[float, np.ndarray]:
    idx = 0 if which == "kappa" else 1

    def neg(v: np.ndarray) -> float:
        x = np.array([[v[0], v[1]], [0.0, v[2]]])
        n = opnorm(x)
        if n < 1e-9:
            return 0.0
        x = x / n
        d = defects(x)[idx]
        if d <= 1e-12:
            return 0.0
        g = D1_exact(x)
        return -(d / g) if g > 1e-12 else -1e9

    v0 = np.array([x0[0, 0], x0[0, 1], x0[1, 1]])
    res = minimize(neg, v0, method="Nelder-Mead",
                   options={"maxiter": 6000, "xatol": 1e-13, "fatol": 1e-15})
    x = np.array([[res.x[0], res.x[1]], [0.0, res.x[2]]])
    return -float(res.fun), x / opnorm(x)


def main() -> None:
    print("validating the dual gauge against the primal ...")
    for x in (np.array([[0.0, 1.0], [0.0, 0.0]]),
              np.diag([1.0, -1.0]),
              np.array([[0.39982, -0.84014], [0.0, -0.39982]])):
        x = x / opnorm(x)
        print(f"  dual {D1_exact(x):.8f}   primal (sampled) {D1_primal(x):.8f}")

    print("\nsweeping T_2(R) ...")
    best = sweep()
    for key in ("kappa", "lambda"):
        val, x = best[key]
        if x is None:
            print(f"  {key} = 0")
            continue
        rval, rx = refine(x, key)
        d = defects(rx)[0 if key == "kappa" else 1]
        print(f"  {key}: grid {val:.8f}   refined {rval:.8f}")
        print(f"      witness =\n{np.round(rx, 6)}")
        print(f"      defect = {d:.8f},  D_1 = {D1_exact(rx):.8f}")


if __name__ == "__main__":
    main()
