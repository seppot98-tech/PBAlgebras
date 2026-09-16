# PBAlgebras

**Perfectly bounded real Banach algebras** — an involution-free extension of
Gelfand duality, and a category built for noncommutative geometry organised by
*normality* rather than by an involution.

---

## The idea in one page

Gelfand duality identifies commutative unital C\*-algebras with compact
Hausdorff spaces. Noncommutative geometry then drops commutativity and keeps the
involution: the C\*-identity `‖a‖² = ‖a*a‖` makes the norm a **singular-value**
radius.

There is another road. On a commutative C\*-algebra one also has

```
‖a²‖ = ‖a‖²                                    (SQ), the square property
```

which makes the norm an **eigenvalue** radius — it says exactly that
`‖a‖ = r(a)`. On a C\*-algebra the two agree precisely on normal elements.
Keeping (SQ) instead of the involution is the fork this project takes.

Taken literally the fork is a dead end: (SQ) together with a spectral-reality
condition **forces commutativity**, so the noncommutative theory would have no
objects. The fix is a *commutator budget* — a quantity `D(a)`, homogeneous of
degree 2, that vanishes exactly on the centre:

```
‖a‖² − ‖a²‖            ≤ C·D(a)                (PB1)  approximate normality
‖ ‖a²‖·1 − a² ‖ − ‖a²‖ ≤ C·D(a)                (PB2)  approximate spectral reality
```

Because the budget is invisible on the centre, the commutative theory is
untouched; off the centre it funds a controlled failure. With `D₁(a) = ‖ad_a‖²`
and `C = ¼` the class contains `C(X,ℝ)`, `ℍ`, `M_n(ℝ)`, `T_n(ℝ)` — and excludes
`ℂ` and `ℝ[ε]/(ε²)`. **Noncommutative points are admitted; infinitesimal ones
are not**, and for the same reason: an unfunded defect.

## Why it is a category worth having

Three structural results (all proved in the paper):

| Result | Statement |
|---|---|
| **The base** | `Z(A) ≅ C(X_A,ℝ)`: every object carries a canonical compact Hausdorff base, *extracted* from its centre, not imposed. `M₂(ℝ)` and `ℍ` are *pointless* — the base is a point — so the theory is not secretly about bundles over a nontrivial space. |
| **Functoriality** | The morphism law forces `f(Z(A)) ⊆ Z(B)` automatically, so `A ↦ X_A` is a functor `PB → CompHaus^op`. This works for `D₁` and **fails for every `D_q`, `q ≥ 2`**. |
| **Coreflection** | `Z` is right adjoint to `CommPB ↪ PB`, and `CommPB ≃ CompHaus^op`. Every PB-algebra is a `C(X_A,ℝ)`-algebra, hence fibred over its own base. |

Plus a sharp threshold: for `C < ¼` every object satisfies `r(a) ≥ (1−4C)‖a‖`,
hence is semisimple with no nilpotents; at `C = ¼` the noncommutative examples
enter all at once. **`¼` is a phase transition, not a normalisation artefact.**

## What is new here, relative to the working notes

- **A claim refuted, with a closed form.** The notes state that no `D₁`-constant
  separates `T₂(ℝ)` from `M₂(ℝ)`: both were said to have critical ratio `¼`,
  attained at the shared element `E₁₂`. False. The extremiser of `T₂(ℝ)` is the
  *traceless* element `[[c, 1−c²],[0,−c]]`, and

  ```
  κ(T₂(ℝ)) ≥ κ* = 4y(1−y)/(1+y)⁴ ,   y = c² = (5−√17)/4 ,   κ* = 0.3098421747…
  ```

  which is `> ¼`. So `T₂(ℝ) ∉ PB^{¼}` while `M₂(ℝ) ∈ PB^{¼}`, and every constant
  in `[¼, κ*)` separates them. Proved (weak duality for the nuclear norm), not
  just computed.

  This matters structurally. The notes present `D₁` vs `D_∞` as a forced trade:
  `D₁` keeps the functorial base but cannot exclude non-central radical, `D_∞`
  excludes it but breaks closure. At `C = ¼` the trade is not forced — which
  raises the **conjecture that every `PB^{¼}` algebra is semisimple**: at
  exactly one constant, the norm gauge would exclude the radical, as the
  C\*-identity does with an involution, while keeping the base functorial.

  Every radical-carrying algebra we can compute is excluded at `¼` —
  `T₂(ℝ)`, `T₃(ℝ)` and both parabolic subalgebras of `M₃(ℝ)`, the last two
  being the informative cases since they carry a full `M₂(ℝ)` block with which
  to fund a radical defect, and still cannot (`make scan`).

- **An open question answered — for large constants.** The notes asked whether a
  legal cone can reach a non-scalar diagonal in `M₂(ℝ)`, which would rule out
  `ℝ1` as the equalizer of `id` and `Ad_diag(1,−1)`. It can: `T₂(ℝ)`, mapped
  onto the diagonal by killing its radical, is such a cone. But a cone must be
  an *object*, and `T₂(ℝ)` is one only for `C ≥ κ*` — so the answer is **yes**
  above `κ*` and **open** at `C = ¼`. The question is then reduced to the
  representability of an explicit functor.

- **Two corrections.** Monomorphisms are **not** the gauge-preserving embeddings
  (any injective morphism is monic; the gauge-preserving embeddings are the
  right notion of *subobject*). And `D_q` does **not** vanish exactly on the
  centre for `q ≥ 2` — `diag(1,−1) ∈ T₂(ℝ)` is a counterexample — which is
  precisely why the base functor needs `D₁`.

- **One programme retired.** Cocompleteness was to be checked via free products.
  But the coproduct of `C(X,ℝ)` and `C(Y,ℝ)` in `PB` is `C(X×Y,ℝ)`, so the free
  product is *not* the coproduct.

## Layout

```
paper/        the mathematical document (LaTeX)
  main.tex      build with `make paper`
  sections/     one file per section
verify/       numerical verification package
  pbalg/        gauges, defects, critical constants
  tests/        the paper's hand computations, as tests
  experiments/  reproducible runs that emit the paper's tables
docs/         status ledger, open problems, terminology map, roadmap
```

## Build

```bash
make paper     # -> paper/main.pdf
make test      # pytest over verify/
make verify    # recompute the constants, regenerate the paper's table
make all
```

Requirements: a TeX installation with `amsmath`/`amsthm`/`mathtools`, and
Python ≥ 3.10 with `numpy`, `scipy`, `pytest`.

## Status discipline

Every substantive claim in the paper carries one of four tags:

- **proved** — a proof is given, or a precise citation supplied;
- **numerical** — a finite-dimensional inequality checked by `verify/`, *not*
  proved. A reported supremum is a lower bound found by search;
- **conjectural** / **open** — what they say.

The ledger is [`docs/STATUS.md`](docs/STATUS.md) and §12 of the paper. The
intent is that the open column shrinks over time and that nothing silently
migrates from *numerical* to *proved*.

A caution that is part of the method. The first search reported
`κ(T₂(ℝ)) = 0.317 > ¼` — a claim that, if true, would falsify a proposition of
the notes. It was *also* partly an artefact: the gauge optimiser was
under-resolving on proper subalgebras, and an underestimated budget inflates
every ratio built on it. Chasing that down produced an exact method for `T₂(ℝ)`
(the gauge there is a support function with a closed dual form, so no search is
involved), which put the true value at `0.3098421747…`. The headline was right
and the first number was wrong, which is the normal situation. It is now a
theorem rather than a number, and the searches refine every witness with a
high-accuracy gauge before reporting it.

## Where to start

- The mathematics: `paper/main.tex`, §§1–2 for the rigid core, §7 for the
  coreflection, §9 for the categorical status.
- The open problems: [`docs/OPEN-PROBLEMS.md`](docs/OPEN-PROBLEMS.md), ordered
  by leverage. Problems 1–3 would change the shape of the theory.
- If you are coming from the working notes: the rename map is in
  [`docs/TERMINOLOGY.md`](docs/TERMINOLOGY.md).
