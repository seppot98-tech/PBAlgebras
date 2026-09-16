"""The identity D_1(a) = 4 dist(a, R1)^2: where it holds and where it fails.

The paper uses it as a fast path for M_n(R) and H, and records that it fails on
triangular algebras.  Both halves are tested here, because the fast path would
silently corrupt the critical-constant search if the identity were false.
"""

import numpy as np
import pytest

from pbalg.algebras import full_matrix, quaternions, upper_triangular
from pbalg.gauges import D1, dist_to_scalars, opnorm

TOL = 1e-6


@pytest.mark.parametrize("alg", [full_matrix(2), full_matrix(3), quaternions()])
def test_identity_holds(alg):
    rng = np.random.default_rng(4)
    for _ in range(8):
        a = alg.element(rng.standard_normal(alg.dim))
        a = a / opnorm(a)
        assert D1(alg, a) == pytest.approx(4 * dist_to_scalars(a) ** 2, abs=1e-5)


@pytest.mark.parametrize("alg", [upper_triangular(2), upper_triangular(3)])
def test_identity_fails_on_triangular(alg):
    """Strictly less, for at least some elements -- and never more."""
    rng = np.random.default_rng(4)
    strict = 0
    for _ in range(8):
        a = alg.element(rng.standard_normal(alg.dim))
        a = a / opnorm(a)
        lhs, rhs = D1(alg, a), 4 * dist_to_scalars(a) ** 2
        assert lhs <= rhs + TOL  # Lemma: D_1(a) <= 4 dist(a, Z(A))^2
        if lhs < rhs - 1e-3:
            strict += 1
    assert strict > 0, "expected the Stampfli identity to fail on T_n(R)"
