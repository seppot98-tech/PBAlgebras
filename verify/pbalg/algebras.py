"""Finite-dimensional real algebras, presented as spans of real matrices.

An :class:`Algebra` is a unital subalgebra of ``M_n(R)`` given by a linear
basis.  The norm is always the operator norm inherited from the Euclidean
structure on ``R^n``; for the quaternions this is realised by the left regular
representation, for which the operator norm is the (multiplicative) quaternion
modulus.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np


@dataclass(frozen=True)
class Algebra:
    """A unital subalgebra of ``M_n(R)`` spanned by ``basis``."""

    name: str
    basis: tuple[np.ndarray, ...]
    latex: str = ""

    @property
    def dim(self) -> int:
        return len(self.basis)

    @property
    def size(self) -> int:
        return self.basis[0].shape[0]

    @property
    def one(self) -> np.ndarray:
        return np.eye(self.size)

    def element(self, coeffs: Sequence[float]) -> np.ndarray:
        """The element with the given coordinates in ``basis``."""
        out = np.zeros_like(self.basis[0])
        for c, b in zip(coeffs, self.basis):
            out = out + c * b
        return out

    def random(self, rng: np.random.Generator) -> np.ndarray:
        return self.element(rng.standard_normal(self.dim))

    def contains(self, m: np.ndarray, tol: float = 1e-9) -> bool:
        """Is ``m`` in the span of the basis?"""
        A = np.stack([b.ravel() for b in self.basis], axis=1)
        sol, *_ = np.linalg.lstsq(A, m.ravel(), rcond=None)
        return bool(np.linalg.norm(A @ sol - m.ravel()) < tol)

    def is_closed_under_multiplication(self, tol: float = 1e-9) -> bool:
        return all(
            self.contains(x @ y, tol) for x in self.basis for y in self.basis
        )


def _unit(n: int, i: int, j: int) -> np.ndarray:
    e = np.zeros((n, n))
    e[i, j] = 1.0
    return e


def full_matrix(n: int) -> Algebra:
    """``M_n(R)``."""
    basis = tuple(_unit(n, i, j) for i in range(n) for j in range(n))
    return Algebra(f"M{n}(R)", basis, latex=rf"M_{{{n}}}(\RR)")


def upper_triangular(n: int) -> Algebra:
    """``T_n(R)``, the upper-triangular matrices."""
    basis = tuple(_unit(n, i, j) for i in range(n) for j in range(i, n))
    return Algebra(f"T{n}(R)", basis, latex=rf"T_{{{n}}}(\RR)")


def diagonal(n: int) -> Algebra:
    """``R^n = C(X, R)`` for ``X`` an ``n``-point space."""
    basis = tuple(_unit(n, i, i) for i in range(n))
    return Algebra(f"R^{n}", basis, latex=rf"\RR^{{{n}}}")


def block_upper_triangular(sizes: Sequence[int]) -> Algebra:
    """Block upper-triangular matrices for the given diagonal block sizes.

    ``block_upper_triangular([1, 1])`` is ``T_2(R)``; ``[1, 2]`` is the
    parabolic subalgebra of ``M_3(R)`` whose radical is non-central *and* whose
    diagonal carries a genuine ``M_2(R)`` block.  These are the sharpest test
    cases for the conjecture that every PB-algebra at ``C = 1/4`` is semisimple:
    the matrix block supplies commutators that the radical of ``T_n(R)`` lacks,
    so if any radical-carrying algebra can pay its way at ``1/4``, one of these
    is a candidate.
    """
    n = sum(sizes)
    starts, off = [], 0
    for s in sizes:
        starts.append(off)
        off += s
    basis = []
    for bi, si in enumerate(sizes):
        for bj, sj in enumerate(sizes):
            if bj < bi:
                continue  # strictly below the block diagonal
            for i in range(starts[bi], starts[bi] + si):
                for j in range(starts[bj], starts[bj] + sj):
                    basis.append(_unit(n, i, j))
    name = "P(" + ",".join(str(s) for s in sizes) + ")"
    return Algebra(name, tuple(basis), latex=rf"P_{{{','.join(map(str, sizes))}}}(\RR)")


def quaternions() -> Algebra:
    """``H``, via the left regular representation on ``R^4``."""
    one = np.eye(4)
    i = np.array(
        [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]], dtype=float
    )
    j = np.array(
        [[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0]], dtype=float
    )
    return Algebra("H", (one, i, j, i @ j), latex=r"\HH")


def complexes() -> Algebra:
    """``C`` regarded as a real algebra, via ``a + bi -> [[a,-b],[b,a]]``."""
    one = np.eye(2)
    i = np.array([[0, -1], [1, 0]], dtype=float)
    return Algebra("C", (one, i), latex=r"\CC")


def dual_numbers() -> Algebra:
    """``R[eps]/(eps^2)``, realised inside ``T_2(R)``."""
    return Algebra(
        "R[eps]/(eps^2)", (np.eye(2), _unit(2, 0, 1)), latex=r"\RR[\varepsilon]/(\varepsilon^2)"
    )


def direct_sum(a: Algebra, b: Algebra, name: str | None = None) -> Algebra:
    """Block-diagonal direct sum, i.e. the product algebra with the sup norm."""
    n, m = a.size, b.size

    def embed_left(x: np.ndarray) -> np.ndarray:
        out = np.zeros((n + m, n + m))
        out[:n, :n] = x
        return out

    def embed_right(y: np.ndarray) -> np.ndarray:
        out = np.zeros((n + m, n + m))
        out[n:, n:] = y
        return out

    basis = tuple(embed_left(x) for x in a.basis) + tuple(
        embed_right(y) for y in b.basis
    )
    return Algebra(name or f"{a.name} (+) {b.name}", basis)


#: The running examples of the paper.
CATALOGUE: dict[str, Callable[[], Algebra]] = {
    "R": lambda: diagonal(1),
    "R^2": lambda: diagonal(2),
    "R^3": lambda: diagonal(3),
    "C": complexes,
    "H": quaternions,
    "M2(R)": lambda: full_matrix(2),
    "M3(R)": lambda: full_matrix(3),
    "T2(R)": lambda: upper_triangular(2),
    "T3(R)": lambda: upper_triangular(3),
    "R[eps]/(eps^2)": dual_numbers,
    "M2(R)(+)M2(R)": lambda: direct_sum(full_matrix(2), full_matrix(2)),
    "H(+)R^2": lambda: direct_sum(quaternions(), diagonal(2)),
    "P(1,2)": lambda: block_upper_triangular([1, 2]),
    "P(2,1)": lambda: block_upper_triangular([2, 1]),
}
