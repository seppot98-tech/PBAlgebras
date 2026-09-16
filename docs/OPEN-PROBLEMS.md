# Open problems

Ordered by leverage. Problems 0–3 would change the shape of the theory; the
rest firm up statements already in use. Section references are to `paper/`.

---

### 0. Semisimplicity at the critical constant (§5)

**Conjecture.** Every `A ∈ PB^{¼}_{D₁}` is semisimple, `rad(A) = 0`.

This is the question opened up by the refutation of the notes' claim that no
`D₁`-constant excludes `T₂(ℝ)`. It does: `κ(T₂(ℝ)) = κ* = 0.30984… > ¼`.

Note it does **not** follow from the phase-transition bound `r(a) ≥ (1−4C)‖a‖`,
which is vacuous at `C = ¼`, and it is not about nilpotent *elements* —
`M₂(ℝ) ∈ PB^{¼}` contains `E₁₂` and is simple.

**Evidence.** Every elementary radical-carrying algebra we can compute is
excluded at `¼`. To exclude one it is enough to exhibit a single element with
ratio above `¼` — no supremum needed — and
`verify/experiments/radical_scan.py` finds them:

| Algebra | witness | ratio |
|---|---|---|
| `ℝ[ε]/(ε²)` | `ε` | `∞` (unfunded) |
| `T₂(ℝ)` | `x_0.45` | `0.30894` |
| `T₃(ℝ)` | `x_0.45` in the `{1,3}` corner | `0.30894` |
| `P₁,₂(ℝ)` | `x_0.45` in the `{1,2}` corner | `0.30190` |
| `P₂,₁(ℝ)` | `x_0.45` in the `{1,3}` corner | `0.30047` |

where `x_c = [[c, 1−c²],[0, −c]]`. The parabolic (block upper-triangular)
algebras are the informative cases: unlike `T_n(ℝ)` they carry a full `M₂(ℝ)`
block on the diagonal, so they have far more test elements available to fund a
radical defect — and still cannot fund it at `¼`.

If true, then at exactly one value of the constant the norm gauge excludes the
radical — what the C\*-identity achieves by means of an involution — while
keeping the functorial classical base that `D_∞` and every `D_q`, `q ≥ 2`,
destroy. The "discrimination versus closure" trade-off dissolves at `C = ¼`.

A counterexample would be nearly as valuable: a `PB^{¼}`-algebra with nonzero
radical would be the first object of the theory that no C\*-algebra can
imitate, and by problem 1 below it is also where to look for a legal cone at the
critical constant.

**First steps.** Compute `C_{D₁}(T_n(ℝ))` for `n ≥ 3` and for the other
elementary radical-carrying algebras; then look for a lower bound on `κ` in
terms of `rad(A)`.

---

### 1. Completeness — does the equalizer exist? (§9)

Decide whether `PB^C_{D₁}` has equalizers. The concrete test case is the pair
`id, Ad_u : M₂(ℝ) ⇉ M₂(ℝ)` with `u = diag(1,−1)`.

**What is settled.** The set-theoretic equalizer is the diagonal `Δ ≅ ℝ²`, whose
inclusion is *not* legal (`D_Δ ≡ 0` while `D_{M₂}(diag(α,β)) = (α−β)²`). The
natural candidate is then the largest saturated subalgebra `ℝ1` — and for
`C ≥ κ*` it is ruled out: `T₂(ℝ)`, mapped onto `Δ` by killing its radical, is a
legal cone with non-scalar image.

**The answer depends on the constant.** A cone must be an object, and `T₂(ℝ)` is
one only for `C ≥ κ* = 0.30984…`. So at `C = ¼` — where `M₂(ℝ)` lives and
`T₂(ℝ)` does not — the question is still **open**. Any cone there must, by the
centre lemma, be an object with two distinct characters over a single point of
its base, with their separation funded by its own commutators. If problem 0
holds, such an object is semisimple, which rules out the obvious construction.

**The sharpened question.** Legal cones into `Δ` correspond to pairs of
characters `(χ₁, χ₂)` with `(χ₁(x) − χ₂(x))² ≤ D₁(x)` for all `x` (call these
*admissible*). So the equalizer exists **iff** the functor

```
F(T) = { admissible pairs of characters of T }
```

is representable. Two constraints a representing object `E` must satisfy:
`F(M₂(ℝ)) = ∅`, so no morphism `M₂(ℝ) → E` exists; and `Hom(C(X,ℝ), E)` must be
`C(X_E, X)`.

A negative answer makes `PB` incomplete — the expected outcome. A positive one
produces the free PB-algebra on a pair of characters whose separation is funded
by its own commutators, which would be an interesting object in its own right.

---

### 2. Closure under quotients (§9, §7.3)

Is `PB^C_{D₁}` closed under quotients by closed two-sided ideals? Quotient maps
are always legal; at issue is whether `A/I` still satisfies the axioms with the
same constant. The risk is that the quotient collapses the test elements
realising a budget faster than it shrinks the defect.

This one question controls two others: **coequalizers**, hence cocompleteness;
and whether the **fibres `A_x`** of the bundle picture are PB-algebras. A
counterexample would be nearly as valuable as a proof — it would name the extra
hypothesis the bundle picture needs.

---

### 3. Coproducts (§9)

Does `PB^C_{D₁}` have coproducts? Note what is *not* the answer: the free
product. The coproduct of `C(X,ℝ)` and `C(Y,ℝ)` is `C(X×Y,ℝ)`, because every
PB-morphism out of a commutative object lands in the centre of its target. A
natural guess is a budget-completed free product. First case to compute:
`M₂(ℝ) ⨿ M₂(ℝ)`.

---

### 4. The uniform C\*-bound (§10)

Prove or refute: for every element of a real C\*-algebra,

```
‖a‖² − ‖a²‖ ≤ dist(a, ℝ1)² .
```

Equivalently (via the Stampfli identity) every real C\*-algebra satisfies (PB1)
at `C = ¼`. Evidence: equality at every square-zero element, triviality at every
normal element, and a numerical search over `M_n(ℝ)`, `n ≤ 4`. A proof would
reduce membership in `𝓘` to the spectral-reality axiom alone — a condition on
the centre.

---

### 5. Real Stampfli identity (§3)

Prove `‖ad_a‖ = 2·dist(a, ℝ1)` for `a ∈ M_n(ℝ)`, or find the correct real
statement. Stampfli's theorem is over `ℂ`. Numerically the identity holds for
`M_n(ℝ)` with `n ≤ 4` and for `ℍ` (proved there), and **fails** on `T_n(ℝ)` —
so whatever is true is sensitive to more than the norm.

---

### 6. Matrix stability (§13)

Is `κ_{D₁}(M_n(ℝ)) = ¼` for all `n`, and is `C_{D₁}(M_n(A))` bounded in `n` for
fixed `A`? Uniform boundedness is the minimum requirement for any K-theoretic or
Morita-theoretic development. Without it, the category cannot carry the usual
noncommutative-geometry machinery. Confirmed numerically for `n ≤ 4`.

---

### 7. Quaternionic Stone–Weierstrass (§5)

Prove that (SQ) together with (PB2) forces an isometric unital embedding into an
algebra of sections with `ℝ`- and `ℍ`-fibres — the infinite-dimensional form of
"the zero-budget survivors are products of `ℝ` and `ℍ`", i.e. Dyson's threefold
way with the unitary class excised.

---

### 8. A finite-`q` sweet spot (§11)

Is there a single finite `q` whose class excludes `T₂(ℝ)` while remaining closed
under quotients? The cost is known in advance: for every `q ≥ 2` the vanishing
locus of `D_q` is strictly larger than the centre, so morphisms no longer
preserve centres and **the base stops being functorial**. Any positive answer
should come with a replacement for §7.

---

### 9. Real Hirschfeld–Żelazko (§2)

Give a self-contained proof that a real Banach algebra with `r = ‖·‖` is
commutative. This is the one citation-dependent step in the rigid-core theorem.
The obstacle: `r = ‖·‖` on `A` does not transfer verbatim to `A_ℂ`.

---

### 10. An intrinsic budget (§9)

Is there an element-intrinsic, involution-free replacement for the gauge — one
*not* defined as a supremum over the ambient algebra — agreeing with `D₁` on the
examples?

This is the most consequential question on the list. The gauge is
ambient-referencing, which is exactly why `PB` is not closed under subalgebras
and why equalizers are hard: an equation can delete the very test elements that
realise a budget. A C\*-algebra has no such failure, because `r(a*a)^{1/2}` is
built from `a` and `a*` alone. An intrinsic budget would be inherited by
subalgebras and the C\*-proof of completeness would go through verbatim. The
answer is plausibly negative — in which case the theorem to prove is that no
such replacement exists.

---

### 11. Faithfulness of the bundle (§7.3)

Is `‖a‖ = sup_x ‖π_x(a)‖`? Equivalently, does the central `C(X_A,ℝ)`-module
structure satisfy the partition-of-unity estimate that holds automatically for
C\*-algebras? A positive answer makes literal the slogan that *a PB-algebra is
the algebra of sections of an upper semicontinuous field of pointless
PB-algebras over its own classical base.*
