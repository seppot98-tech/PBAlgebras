# Status ledger

Mirrors §12 of the paper. **proved** = proof given or precise citation;
**numerical** = checked by `verify/`, not proved; **conjectural** / **open** =
what they say. Entries marked *withdrawn*, *refuted* or *qualified* are claims
of the working notes that did not survive.

## Proved

| Statement | Where |
|---|---|
| (SQ)+(SR) ⟹ `A ≅ C(X,ℝ)`; commutativity is forced, not assumed | §2 † |
| `Z(A) ≅ C(X_A,ℝ)` for every PB-algebra — the classical base | §5 |
| Characters are approximately real: `(Im χ(a))² ≤ C·D₁(a)` | §5 |
| `r(a) ≥ (1−4C)‖a‖`; subcritical objects are semisimple, nilpotent-free | §5 |
| Finite-dimensional zero-budget survivors are products of `ℝ` and `ℍ` | §5 |
| `κ(T₂(ℝ)) ≥ κ* = 0.3098421747… > ¼`, so a constant separates `T₂` from `M₂` | §8 |
| `D₁` vanishes exactly on the centre; `D_q` (`q ≥ 2`) does not | §3 |
| Isometric PB-morphisms preserve the gauge | §6 |
| PB-morphisms preserve centres; `A ↦ X_A` is a functor | §7 |
| `CommPB ≃ CompHaus^op` is **coreflective** in `PB` | §7 |
| Coproducts of commutative objects are `C(X×Y,ℝ)`; the free product is not the coproduct | §7 |
| `PB^C` has all products; the union over `C` does not | §9 |
| No terminal object without the degenerate algebra | §9 |
| `PB` is not closed under closed subalgebras | §9 |
| `Δ ↪ M₂(ℝ)` is not legal; `ℝ1` is the largest saturated subalgebra of `Δ` | §9 |
| A legal cone onto all of `Δ` exists for `C ≥ κ*`, so `ℝ1` is not the equalizer there | §9 |
| `T₂(ℝ)` has non-central radical, lies outside `𝓘` | §8, §11 |
| `D_∞` excludes `T₂(ℝ)` and admits `M₂(ℝ)` | §11 |
| `U` is faithful, bijective on objects, not full | §10 |
| Chirally asymmetric C\*-algebras have no real form | §8 |

† modulo a real-scalar Hirschfeld–Żelazko input; see problem 9.

## Numerical

| Statement | Where |
|---|---|
| `D₁(a) = 4·dist(a,ℝ1)²` on `M_n(ℝ)`, `n ≤ 4` (proved for `ℍ`) | §3 |
| The same identity **fails** on `T_n(ℝ)` | §3 |
| `κ = λ = ¼` for `M₂(ℝ)`; `κ = 0`, `λ = ¼` for `ℍ` | §8 |
| `κ(T₂(ℝ)) = κ*` exactly, i.e. the traceless slice is extremal | §8 |
| `κ(M_n(ℝ)) = ¼` for `n ≤ 4` | §13 |

## Conjectural

| Statement | Where |
|---|---|
| Every `PB^{¼}_{D₁}`-algebra is semisimple | §5, problem 0 |
| Uniform C\*-bound: `‖a‖² − ‖a²‖ ≤ dist(a,ℝ1)²` on real C\*-algebras | §10, problem 4 |
| Quaternionic Stone–Weierstrass | §5, problem 7 |

## Open

| Statement | Where |
|---|---|
| Completeness / representability of the admissible-pair functor | problem 1 |
| The cone question at `C ∈ [¼, κ*)`, in particular at the critical constant | §9 |
| Closure under quotients (⟹ coequalizers, and the bundle fibres) | problem 2 |
| Coproducts, hence cocompleteness | problem 3 |
| Faithfulness of the bundle | problem 11 |
| Existence of an intrinsic budget | problem 10 |
| A finite-`q` sweet spot | problem 8 |

## Did not survive

| Claim of the working notes | Fate |
|---|---|
| No `D₁`-constant separates `T₂(ℝ)` from `M₂(ℝ)`; both have `κ = ¼` at `E₁₂` | **refuted** — `κ(T₂) = κ* > ¼`, extremiser traceless (§8) |
| Monomorphisms are the gauge-preserving embeddings | **withdrawn** — any injective morphism is monic; these are the *subobjects* (§6) |
| Every `D_q` vanishes exactly on the centre | **withdrawn** — `diag(1,−1) ∈ T₂(ℝ)` (§3) |
| Cocompleteness via free products | **retired** — the free product is not the coproduct (§7) |
| The excess of `PB` over `𝓘` is exactly non-central radical | **qualified** — false at `C = ¼` as stated (§11) |
| `W` is incomplete because `Δ ↪ M₂` is illegal | already retracted in the notes; the retraction stands (§9) |
| A reflector `BanAlg → PB` gives cocompleteness | already retracted in the notes; it would force completeness (§9) |
