"""Tests of statements the paper proves, checked on the finite-dimensional
examples.  These are regression tests for the *theory*, not for the code: each
one would fail if the corresponding proposition were false.
"""

import numpy as np
import pytest

from pbalg.algebras import (
    diagonal,
    direct_sum,
    dual_numbers,
    full_matrix,
    quaternions,
    upper_triangular,
)
from pbalg.defects import pb1_defect, pb2_defect
from pbalg.gauges import D1, opnorm
from pbalg.search import critical_constants

TOL = 1e-5


def test_products_take_the_sup_of_the_gauge():
    """Proposition "Products": D(a) = sup_i D(a_i)."""
    a1 = np.array([[0.0, 1.0], [0.0, 0.0]])          # D_1 = 4 in M_2
    a2 = np.diag([1.0, -1.0])                        # D_1 = 4 in M_2
    prod = direct_sum(full_matrix(2), full_matrix(2))
    for x, y, expected in ((a1, 0 * a1, 4.0), (0 * a1, a2, 4.0), (a1, a2, 4.0)):
        block = np.zeros((4, 4))
        block[:2, :2] = x
        block[2:, 2:] = y
        assert D1(prod, block) == pytest.approx(expected, abs=1e-4)


def test_excluded_algebras_are_excluded():
    """Dual numbers and C have unfunded defects: not PB-algebras at any C."""
    eps = np.array([[0.0, 1.0], [0.0, 0.0]])
    assert D1(dual_numbers(), eps) == pytest.approx(0.0, abs=TOL)
    assert pb1_defect(eps) > 0.5              # unfunded square defect

    i = np.array([[0.0, -1.0], [1.0, 0.0]])
    from pbalg.algebras import complexes
    assert D1(complexes(), i) == pytest.approx(0.0, abs=TOL)
    assert pb2_defect(i) > 0.5                # unfunded reality defect


def test_commutative_objects_have_no_defects():
    rng = np.random.default_rng(2)
    alg = diagonal(4)
    for _ in range(10):
        a = alg.element(rng.standard_normal(alg.dim))
        assert pb1_defect(a) <= TOL
        assert pb2_defect(a) <= TOL


def test_subalgebras_need_not_be_pb_algebras():
    """R[eps]/(eps^2) is a closed unital subalgebra of M_2(R)."""
    alg = dual_numbers()
    assert alg.is_closed_under_multiplication()
    r = critical_constants(alg, tries=4)
    assert not r.is_pb


def test_the_T2_cone_is_a_legal_morphism():
    """The cone of the equalizer theorem: t(x) = diag(chi_1 x, chi_2 x).

    Legality is the inequality (chi_1 - chi_2)^2 <= D_1^{T_2}, which the proof
    establishes with the witness z = E_12.
    """
    T2, M2 = upper_triangular(2), full_matrix(2)
    rng = np.random.default_rng(3)
    for _ in range(25):
        x = T2.element(rng.standard_normal(T2.dim))
        alpha, beta = x[0, 0], x[1, 1]
        t_x = np.diag([alpha, beta])
        assert opnorm(t_x) <= opnorm(x) + TOL                  # contractive
        assert D1(M2, t_x) <= D1(T2, x) + TOL                  # budget-contracting
        assert D1(M2, t_x) == pytest.approx((alpha - beta) ** 2, abs=1e-5)


def test_M2_and_H_have_critical_constant_one_quarter():
    for alg in (full_matrix(2), quaternions()):
        r = critical_constants(alg, tries=12)
        assert r.constant == pytest.approx(0.25, abs=2e-3), r.summary()
