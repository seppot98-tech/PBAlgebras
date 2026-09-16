"""The hand computations of the paper, as tests.

Every number here is one that the text asserts; if a change to the gauge code
breaks one of these, the text is wrong or the code is.
"""

import numpy as np
import pytest

from pbalg.algebras import (
    complexes,
    diagonal,
    dual_numbers,
    full_matrix,
    quaternions,
    upper_triangular,
)
from pbalg.defects import pb1_defect, pb2_defect
from pbalg.gauges import D1, Dinf, Dq, dist_to_scalars, opnorm

M2 = full_matrix(2)
T2 = upper_triangular(2)
H = quaternions()

E12 = np.array([[0.0, 1.0], [0.0, 0.0]])
U = np.diag([1.0, -1.0])
J = np.array([[0.0, -1.0], [1.0, 0.0]])

TOL = 1e-6


def test_catalogue_algebras_are_algebras():
    for alg in (M2, T2, H, diagonal(3), dual_numbers(), complexes(),
                full_matrix(3), upper_triangular(3)):
        assert alg.is_closed_under_multiplication(), alg.name


def test_E12_gauges():
    # Example "The fundamental computation" in the paper.
    assert D1(M2, E12) == pytest.approx(4.0, abs=TOL)
    assert D1(T2, E12) == pytest.approx(4.0, abs=TOL)
    assert Dinf(M2, E12) == pytest.approx(1.0, abs=TOL)
    assert Dinf(T2, E12) == pytest.approx(0.0, abs=TOL)
    assert Dq(M2, E12, 2) == pytest.approx(1.0, abs=TOL)


def test_E12_is_the_extremiser_at_one_quarter():
    # defect 1, budget 4: the ratio 1/4 of the phase transition.
    assert pb1_defect(E12) == pytest.approx(1.0, abs=TOL)
    assert pb1_defect(E12) / D1(M2, E12) == pytest.approx(0.25, abs=TOL)


def test_rotation_realises_the_pb2_defect():
    # J^2 = -1: reality defect 1, budget 4.
    assert pb2_defect(J) == pytest.approx(1.0, abs=TOL)
    assert D1(M2, J) == pytest.approx(4.0, abs=TOL)
    assert pb2_defect(J) / D1(M2, J) == pytest.approx(0.25, abs=TOL)


def test_quaternion_unit_realises_the_pb2_defect():
    i = H.basis[1]
    assert pb1_defect(i) == pytest.approx(0.0, abs=TOL)  # norm is multiplicative
    assert pb2_defect(i) == pytest.approx(1.0, abs=TOL)
    assert D1(H, i) == pytest.approx(4.0, abs=TOL)


def test_diagonal_gauge_in_M2():
    # D_1(diag(a,b)) = (a-b)^2, used in the equalizer discussion.
    for a, b in ((1.0, -1.0), (2.0, 0.5), (0.3, 0.3)):
        d = np.diag([a, b])
        assert D1(M2, d) == pytest.approx((a - b) ** 2, abs=TOL)


def test_gauge_ordering_and_universal_bound():
    rng = np.random.default_rng(0)
    for _ in range(5):
        a = M2.element(rng.standard_normal(M2.dim))
        a = a / opnorm(a)
        d1, d2, dinf = D1(M2, a), Dq(M2, a, 2), Dinf(M2, a)
        assert dinf <= d2 + TOL
        assert d2 <= d1 + TOL
        assert d1 <= 4.0 + TOL


def test_gauge_is_translation_invariant_and_homogeneous():
    rng = np.random.default_rng(1)
    a = M2.element(rng.standard_normal(M2.dim))
    base = D1(M2, a)
    assert D1(M2, a + 3.0 * np.eye(2)) == pytest.approx(base, rel=1e-6)
    assert D1(M2, 2.0 * a) == pytest.approx(4.0 * base, rel=1e-6)


def test_gauge_vanishes_exactly_on_the_centre_for_D1():
    # central elements: zero budget
    assert D1(M2, np.eye(2)) == pytest.approx(0.0, abs=TOL)
    # but D_q for q >= 2 vanishes off the centre in T_2 -- Warning in the paper
    assert D1(T2, U) == pytest.approx(4.0, abs=TOL)
    assert Dq(T2, U, 2) == pytest.approx(0.0, abs=TOL)
    assert Dinf(T2, U) == pytest.approx(0.0, abs=TOL)


def test_dist_to_scalars():
    assert dist_to_scalars(E12) == pytest.approx(1.0, abs=TOL)
    assert dist_to_scalars(np.eye(2)) == pytest.approx(0.0, abs=TOL)
    assert dist_to_scalars(np.diag([1.0, -1.0])) == pytest.approx(1.0, abs=TOL)
