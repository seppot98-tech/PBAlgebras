# Terminology and notation

## Rename map from the working notes

| Working notes | Here | Why |
|---|---|---|
| `W`-algebra, `W`, `W_q`, `W^C_q` | **PB-algebra**, `PB`, `PB_D`, `PB^C_D` | "W-algebra" collides with W\*-algebra (von Neumann algebra); the whole point is that there is no involution here, so the collision is actively misleading. "Perfectly bounded" is the author's name for the objects and is adopted as given. A reading consistent with the axioms, offered as a gloss rather than a definition: the norm is *perfectly* pinned to the spectral radius where the budget vanishes, and *boundedly* away from it elsewhere. |
| `(C2)` | **(SR)**, spectral reality | "Condition 2" names nothing. On a commutative algebra the condition says exactly that no character takes a non-real value. |
| `W1`, `W2` | **(PB1)**, **(PB2)** | Consistency with the object name. |
| `W1₀` | **`PB⁰`** | The unrelaxed (zero-budget) case. |
| `†`-morphism | **PB-morphism** | `†` is standard for dagger categories, where it is an involutive contravariant identity-on-objects functor. This notion is neither involutive nor contravariant, and the theory has no involution. |
| `D`-saturated subalgebra | **saturated subalgebra**; the inclusion is a **saturated embedding** | Shorter, and the gauge is clear from context. |
| classical base | **classical base** `X_A = Spec Z(A)` (unchanged) | |
| pointless object | **pointless object** (unchanged) | `Z(A) = ℝ1`, i.e. the base is a single point. |
| the intersection `I` | **`𝓘`** (unchanged in meaning) | Real forms of C\*-algebras satisfying (PB1), (PB2). |

## Notation

| Symbol | Meaning |
|---|---|
| `A` | a unital real Banach algebra, `‖1‖ = 1`, norm submultiplicative |
| `A_ℂ` | complexification |
| `Z(A)`, `rad(A)` | centre, Jacobson radical |
| `r(a)` | spectral radius |
| `ad_a` | `[a, −]` |
| `D₁(a)` | `‖ad_a‖²`, the norm gauge |
| `D_q(a)` | `sup_{‖z‖≤1} ‖[a,z]^q‖^{2/q}` |
| `D_∞(a)` | `sup_{‖z‖≤1} r([a,z])²`, the spectral gauge |
| `κ_D(A)`, `λ_D(A)` | critical constants for (PB1), (PB2) |
| `C_D(A)` | `max(κ, λ)`; `A ∈ PB^C` iff `C_D(A) ≤ C` |
| `X_A` | classical base |
| `PB^C_D` | the category, for gauge `D` and constant `C` |
| `CommPB` | full subcategory of commutative objects, `≃ CompHaus^op` |

## Conventions worth knowing

- **The degenerate object.** The zero algebra (`1 = 0`) is admitted. Without it
  there is no terminal object and hence no limits, for reasons having nothing to
  do with the subject. Unital C\*-algebra theory makes the same convention.
- **Gauges are translation invariant**: `D(a + λ1) = D(a)`, since `ad_{λ1} = 0`.
  The *defects* are not, which is why several proofs translate by a scalar
  before applying an axiom.
- **The constant is a filtration, not an invariant.** `C_D(A)` is the invariant;
  `PB^C` collects the objects below a threshold. Products need a *uniform* `C`,
  which is why `PB^C`, and not the union over `C`, is the category of interest.
