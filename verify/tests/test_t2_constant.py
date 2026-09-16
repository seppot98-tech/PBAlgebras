"""The critical constant of T_2(R), which refutes a claim of the working notes.

The paper proves kappa(T_2(R)) >= kappa* = 4y(1-y)/(1+y)^4 at y = (5-sqrt17)/4,
via a rank-one nuclear-norm dual certificate.  These tests pin every step of
that proof, and the conclusion kappa* > 1/4.
"""

import math

import numpy as np
import pytest

from experiments.t2_exact import D1_exact, D1_primal, defects
from pbalg.algebras import upper_triangular
from pbalg.gauges import D1, opnorm

Y_STAR = (5 - math.sqrt(17)) / 4
C_STAR = math.sqrt(Y_STAR)
KAPPA_STAR = 4 * Y_STAR * (1 - Y_STAR) / (1 + Y_STAR) ** 4

T2 = upper_triangular(2)


def x_c(c: float) -> np.ndarray:
    return np.array([[c, 1 - c * c], [0.0, -c]])


def test_kappa_star_value_and_closed_form():
    assert KAPPA_STAR == pytest.approx(0.3098421747, abs=1e-9)
    closed = 128 * (3 * math.sqrt(17) - 11) / (15112 - 3528 * math.sqrt(17))
    assert closed == pytest.approx(KAPPA_STAR, rel=1e-12)


def test_kappa_star_exceeds_one_quarter():
    """The refutation: a constant in [1/4, kappa*) separates T_2 from M_2."""
    assert KAPPA_STAR > 0.25


@pytest.mark.parametrize("c", [0.1, 0.3, C_STAR, 0.7, 0.95])
def test_traceless_family_is_normalised(c):
    """x_c = [[c, 1-c^2],[0,-c]] has norm 1, and x_c^2 = c^2 * 1."""
    x = x_c(c)
    assert opnorm(x) == pytest.approx(1.0, abs=1e-12)
    assert np.allclose(x @ x, c * c * np.eye(2), atol=1e-12)
    assert defects(x)[0] == pytest.approx(1 - c * c, abs=1e-12)


@pytest.mark.parametrize("c", [0.45, C_STAR, 0.6, 0.8])
def test_dual_certificate_bound(c):
    """D_1(x_c) <= (1+c^2)^4/(4c^2), with equality above c = sqrt(2)-1."""
    x = x_c(c)
    bound = (1 + c * c) ** 4 / (4 * c * c)
    # The bound is exact; D1_exact minimises over g numerically, so it can
    # land a few ulps either side of the analytic value.
    assert D1_exact(x) <= bound * (1 + 1e-6)
    assert D1_exact(x) == pytest.approx(bound, rel=1e-6)  # attained here


def test_bound_is_not_attained_below_the_crossover():
    """Below c = sqrt(2)-1 the rank-one certificate is not optimal."""
    c = 0.2
    x = x_c(c)
    assert D1_exact(x) == pytest.approx(4 * (1 - c * c) ** 2, rel=1e-9)
    assert D1_exact(x) < (1 + c * c) ** 4 / (4 * c * c) - 1e-3


def test_optimum_realises_kappa_star():
    x = x_c(C_STAR)
    ratio = defects(x)[0] / D1_exact(x)
    assert ratio == pytest.approx(KAPPA_STAR, rel=1e-7)


def test_dual_agrees_with_primal_and_with_the_general_search():
    """Three independent routes to the same gauge.

    D1_primal is a crude random sampler, so it is a lower bound that converges
    slowly: we require only that it does not exceed the dual value and gets
    within a few percent.  The general search is a real optimiser and must
    agree closely.
    """
    for c in (0.3, C_STAR, 0.7):
        x = x_c(c)
        dual = D1_exact(x)
        sampled = D1_primal(x, samples=60000)
        assert sampled <= dual * (1 + 1e-9)
        assert sampled == pytest.approx(dual, rel=3e-2)
        assert D1(T2, x, restarts=16, iters=80, polish=8) == pytest.approx(
            dual, rel=1e-5
        )


def test_E12_is_not_the_extremiser_of_T2():
    """The error in the working notes: E_12 gives only 1/4, not the supremum."""
    e12 = np.array([[0.0, 1.0], [0.0, 0.0]])
    assert defects(e12)[0] / D1_exact(e12) == pytest.approx(0.25, abs=1e-9)
    assert defects(x_c(C_STAR))[0] / D1_exact(x_c(C_STAR)) > 0.25 + 1e-3
