"""Evidence for the conjecture that every PB-algebra at C = 1/4 is semisimple.

A full critical-constant search is expensive: the outer maximisation over
elements calls the gauge, which is itself a maximisation.  But to *exclude* an
algebra from ``PB^{1/4}`` we do not need the supremum -- one element with

    (||a||^2 - ||a^2||) / D_1(a)  >  1/4

settles it.  So this script hunts for a single witness in each
radical-carrying algebra, which costs one accurate gauge evaluation per
candidate rather than a nested search.

Candidates are the traceless family that is extremal in ``T_2(R)``, embedded in
the larger algebra in every available way, together with random elements and a
short local climb from the best of them.

A witness found is a proof sketch, not a proof: the gauge is computed by search,
so a witness should be re-checked (and, as for ``T_2(R)``, ideally replaced by a
dual certificate) before it is quoted as a theorem.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pbalg.algebras import (  # noqa: E402
    Algebra,
    block_upper_triangular,
    upper_triangular,
)
from pbalg.defects import pb1_defect, pb2_defect  # noqa: E402
from pbalg.gauges import D1, _from_coeffs, _orthonormal_basis, opnorm  # noqa: E402

THRESHOLD = 0.25


def gauge_coarse(alg: Algebra, a: np.ndarray, rng) -> float:
    """Cheap gauge: no polish stage.  It can only *under*estimate ``D_1``."""
    return D1(alg, a, restarts=6, iters=60, polish=0, rng=rng)


def gauge_accurate(alg: Algebra, a: np.ndarray, rng) -> float:
    """High-effort gauge, used only on the handful of screened survivors.

    The effort matters.  At ``polish=4`` this reported a ratio of 0.3064 at a
    P(1,2) candidate whose settled ratio is 0.2939 -- an underestimated budget
    inflates the ratio, the same failure mode that produced a spurious
    ``kappa(T_2(R)) = 0.317`` earlier in this project.  (The scan's reported
    best for P(1,2), 0.3019, is attained at a different element; both exceed
    1/4, so the verdict was never in doubt, but the numbers were.)  Values at
    this setting are stable to ~1e-4 under a further tripling of the effort.
    """
    return D1(alg, a, restarts=30, iters=120, polish=14, rng=rng)


def ratios(alg: Algebra, a: np.ndarray, rng, accurate: bool) -> tuple[float, float]:
    """``(kappa-ratio, lambda-ratio)`` at ``a``.

    With ``accurate=False`` the budget is underestimated, so the ratios are
    *over*estimates.  That is the right direction for screening: an element
    whose coarse ratio is below the threshold cannot be a witness.
    """
    n = opnorm(a)
    if n < 1e-12:
        return 0.0, 0.0
    a = a / n
    g = (gauge_accurate if accurate else gauge_coarse)(alg, a, rng)
    d1, d2 = pb1_defect(a), max(pb2_defect(a), 0.0)
    if g <= 1e-9:
        return (float("inf") if d1 > 1e-9 else 0.0,
                float("inf") if d2 > 1e-9 else 0.0)
    return d1 / g, d2 / g


def traceless_embeddings(alg: Algebra, c: float) -> list[np.ndarray]:
    """The T_2 extremiser [[c, 1-c^2],[0,-c]] placed in every 2x2 corner."""
    n = alg.size
    out = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            m = np.zeros((n, n))
            m[i, i] = c
            m[j, j] = -c
            m[i, j] = 1 - c * c
            if alg.contains(m):
                out.append(m)
    return out


def scan(alg: Algebra, seed: int = 0, randoms: int = 80,
         climb_iter: int = 250, verify_top: int = 12) -> dict:
    """Screen cheaply, then verify the best candidates with the accurate gauge."""
    rng = np.random.default_rng(seed)
    q = _orthonormal_basis(alg)
    n = alg.size

    candidates: list[np.ndarray] = []
    for c in np.linspace(0.2, 0.8, 13):
        candidates.extend(traceless_embeddings(alg, float(c)))
    candidates.extend(alg.element(rng.standard_normal(alg.dim))
                      for _ in range(randoms))

    # a local climb on the cheap objective, started from the best screened point
    screened = [(ratios(alg, a, rng, accurate=False)[0], a) for a in candidates]
    screened.sort(key=lambda t: -t[0] if np.isfinite(t[0]) else -1e18)
    if screened and screened[0][1] is not None:
        def neg(coef: np.ndarray) -> float:
            return -ratios(alg, _from_coeffs(q, coef, n), rng, accurate=False)[0]

        res = minimize(neg, q.T @ screened[0][1].ravel(), method="Nelder-Mead",
                       options={"maxiter": climb_iter, "xatol": 1e-7,
                                "fatol": 1e-9})
        climbed = _from_coeffs(q, res.x, n)
        if opnorm(climbed) > 1e-12:
            screened.insert(0, (-float(res.fun), climbed))

    best = {"kappa": (0.0, None), "lambda": (0.0, None)}
    for _, a in screened[:verify_top]:
        k, l = ratios(alg, a, rng, accurate=True)
        if k > best["kappa"][0]:
            best["kappa"] = (k, a / opnorm(a))
        if l > best["lambda"][0]:
            best["lambda"] = (l, a / opnorm(a))
    return best


def main() -> None:
    names = sys.argv[1:]
    catalogue = {
        "T2": upper_triangular(2),
        "T3": upper_triangular(3),
        "P12": block_upper_triangular([1, 2]),
        "P21": block_upper_triangular([2, 1]),
    }
    algebras = ([catalogue[n] for n in names] if names
                else list(catalogue.values()))
    print(f"looking for elements with ratio > {THRESHOLD} "
          f"(a witness excludes the algebra from PB^(1/4))\n", flush=True)
    for alg in algebras:
        t0 = time.time()
        best = scan(alg)
        k, kw = best["kappa"]
        l, _ = best["lambda"]
        verdict = "EXCLUDED at C=1/4" if max(k, l) > THRESHOLD else "no witness found"
        print(f"{alg.name:10s} dim={alg.dim:<2d} "
              f"best kappa-ratio={k:.6f}  best lambda-ratio={l:.6f}   "
              f"{verdict}   [{time.time()-t0:.1f}s]", flush=True)
        if kw is not None and k > THRESHOLD:
            print("   witness:")
            for row in np.round(kw / opnorm(kw), 4):
                print("     ", row)


if __name__ == "__main__":
    main()
