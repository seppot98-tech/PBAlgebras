"""Numerical verification tools for the PBAlgebras project.

The package computes, for finite-dimensional real algebras given by a linear
basis of real matrices:

* the commutator gauges ``D_1``, ``D_q`` and ``D_inf``;
* the defects appearing in axioms (PB1) and (PB2);
* the critical constants ``kappa`` and ``lambda`` of an algebra.

Nothing here proves anything.  The point is to turn statements that earlier
drafts carried as "heuristic" into statements backed by a reproducible
optimisation, and to catch false conjectures early.
"""

from .algebras import (
    Algebra,
    full_matrix,
    upper_triangular,
    quaternions,
    diagonal,
    dual_numbers,
    complexes,
    direct_sum,
    CATALOGUE,
)
from .gauges import D1, Dq, Dinf, dist_to_scalars, opnorm
from .defects import pb1_defect, pb2_defect, pb1_ratio, pb2_ratio
from .search import critical_constants

__all__ = [
    "Algebra",
    "full_matrix",
    "upper_triangular",
    "quaternions",
    "diagonal",
    "dual_numbers",
    "complexes",
    "direct_sum",
    "CATALOGUE",
    "D1",
    "Dq",
    "Dinf",
    "dist_to_scalars",
    "opnorm",
    "pb1_defect",
    "pb2_defect",
    "pb1_ratio",
    "pb2_ratio",
    "critical_constants",
]
