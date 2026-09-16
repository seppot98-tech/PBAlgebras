# Roadmap

What it would take to turn this from a well-posed proposal into a theory, in
rough dependency order. Each item says what it unblocks.

## Phase 1 — make the category safe to work in

The current state: `PB^C` has products, an initial and a terminal object, a
functorial classical base and a coreflective commutative part. What is missing
is everything to do with *quotients and subobjects*.

1. **Quotient closure** (problem 2). Until this is settled, the category has no
   coequalizers, the bundle fibres are not known to be objects, and no
   construction that quotients is available. This is the single cheapest
   blocker to remove, and a counterexample would be as useful as a proof.
2. **Semisimplicity at `C = ¼`** (problem 0). Decides whether the theory at its
   critical constant is a radical-free theory. Shapes everything downstream: if
   true, the objects are much closer to real C\*-algebras than the notes
   expected, and the interesting objects are the *non*-C\* semisimple ones.
3. **Completeness** (problem 1). Expected negative. Worth settling because a
   negative answer redirects effort toward the colimit side and toward problem
   10 (an intrinsic budget), which is the only route to a limit theory.

## Phase 2 — make it a geometry

4. **The bundle picture, faithfully** (problem 11, needs 1). Gives the slogan
   its literal meaning: a PB-algebra is the section algebra of an upper
   semicontinuous field of pointless PB-algebras over its own base. This is the
   structural statement that would make "noncommutative geometry over
   `CompHaus`" more than an analogy.
5. **Classify the pointless objects.** What are the PB-algebras with `Z(A) = ℝ1`?
   In finite dimensions the zero-budget ones are `ℝ` and `ℍ`
   (Dyson, minus the unitary class); at `C = ¼` we know `M_n(ℝ)` and `ℍ` are in.
   A classification of the fibres is what the bundle picture needs to be useful.
6. **Matrix stability** (problem 6). `M_n(A)` must stay in `PB^C` with `C`
   independent of `n`, or no K-theory is possible.

## Phase 3 — invariants

7. **A K-theory.** Requires 6, and an appropriate notion of projection or
   idempotent. Note the theory has no involution, so there is no "self-adjoint
   idempotent"; the natural substitute is an idempotent with a bound on
   `dist(e, Z(A))`, which is a budget condition. Worth checking whether the
   Grothendieck group of such idempotents is homotopy invariant.
8. **Morita theory.** What is the right notion of bimodule when morphisms must
   contract a budget? A PB-bimodule presumably needs its own gauge.
9. **Relation to the tenfold way.** The zero-budget finite-dimensional theory is
   the orthogonal ⊕ symplectic part of Dyson's threefold way. Is the `C = ¼`
   theory the natural home for real-structure-carrying observable algebras? This
   is the physical motivation of §8 and deserves a precise statement.

## Phase 4 — the hard question

10. **An intrinsic budget** (problem 10). Everything awkward about the category
    traces to the gauge being a supremum over the ambient algebra. Either find
    an intrinsic replacement, or prove none exists. A proof of impossibility
    would be a genuine theorem about why this theory cannot have limits, and it
    would explain the contrast with C\*-algebras structurally rather than by
    example.

## Things that are cheap and worth doing early

- Compute `C_{D₁}` for more finite-dimensional algebras, especially
  radical-carrying ones (feeds problem 0 directly). `verify/` already does this;
  it is a matter of adding algebras to `pbalg/algebras.py`.
- Prove the real Stampfli identity (problem 5), which would make the fast path
  in `verify/` a theorem rather than a numerically-checked assumption.
- Close the real Hirschfeld–Żelazko gap (problem 9), the one citation-dependent
  step in the rigid-core theorem.
