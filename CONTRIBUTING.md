# Contributing

This is a mathematics repository with a verification package attached. Two
rules matter more than the rest.

## 1. Tag every claim

Each substantive statement carries exactly one of:

- **proved** — a proof is in the paper, or a precise citation is supplied. If
  the proof leans on an uncited or unverified input, say so in a remark and add
  a problem to `docs/OPEN-PROBLEMS.md`;
- **numerical** — checked by `verify/`, not proved. Say which script and at what
  tolerance;
- **conjectural** / **open**.

Nothing migrates from *numerical* to *proved* without a proof, and the ledger in
`docs/STATUS.md` and §12 of the paper must agree with the text.

## 2. Numerics are searches, not oracles

Every constant in this project is a supremum of a nonsmooth, nonconvex ratio, so
a computed value is a **lower bound found by search**. Two failure modes have
already bitten:

- the *gauge itself* is computed by maximisation, and an under-resolved gauge
  inflates every ratio built on it. An early run reported
  `κ(T₂(ℝ)) = 0.317 > ¼` for exactly this reason, which would have falsified a
  proposition. Always refine a witness with a high-accuracy gauge before
  believing it;
- a fast path can be false. `D₁(a) = 4·dist(a,ℝ1)²` holds on `M_n(ℝ)` and `ℍ`
  but **fails** on `T_n(ℝ)`; `verify/tests/test_stampfli.py` pins both halves.

Where an exact method exists, use it. For `T₂(ℝ)` the gauge is a support
function and has a closed dual form (`experiments/t2_exact.py`) — no search at
all.

## Adding to the paper

- One file per section in `paper/sections/`, included from `main.tex`.
- New results go in the ledger (§12) in the same commit.
- `make paper` must build with no unresolved references; CI enforces this.

## Adding to `verify/`

- New algebras go in `pbalg/algebras.py` and must pass
  `is_closed_under_multiplication`.
- Any hand computation quoted in the paper should appear as a test in
  `verify/tests/`, so that the text and the code cannot drift apart.
- `make test` before pushing.
