# Background page — content outline

A living outline for `pages/3_Background.py`. The page is a crash course on the
math underlying the project, built as a sequence of top-level expander sections,
each containing tabs that pair a short conceptual write-up with an interactive
applet.

## Narrative arc

The page tells four stories, in this order:

1. **The lattice pictures.** What it means to "draw an elliptic curve over
   𝔽ₚ using a lattice": review the basics, then explain that we lift Frobenius
   to a lattice model (Deuring) and draw the picture from that. We *take the
   existence of the lattice for granted* here and defer "how we find it" to
   Story 4.
2. **Isogenies.** Compare the analytic and algebraic descriptions; the analytic
   picture (folding a fundamental domain) is far more intuitive. Demonstrate
   Vélu. Then isogenies over 𝔽ₚ and the volcano structure.
3. **Modular tools** and how we use them to compute isogenies.
4. **The CM bijection** — how we actually find the lattice that Story 1 took for
   granted.

This deliberately shows the *payoff* (pictures) before the *machinery that
justifies it* (the CM bijection). That ordering is intentional.

## Sections

Status legend: ✅ done · ⬜ todo · ↪ moved from elsewhere

The sections are grouped on the page under two `st.header` banners:
**"Visualizing Elliptic Curves"** (§1–§3) and **"Equivalence of Categories"**
(§4 Isogenies onward; future §5–§6 live here too). Miscellaneous topics sits
last.

## Visualizing Elliptic Curves

### §1 — Elliptic curves: the basics  *(foundation / review)*
- ✅ Algebraic curves over ℝ
- ✅ Real elliptic curves + group law
- ✅ Elliptic curves over 𝔽ₚ + group law
- ✅↪ Elliptic curves over ℂ: complex tori *(moved here from "Analytic Methods")*

### §2 — CM Lattices and Frobenius  *(Story 1)*
- ✅↪ Endomorphisms & Complex Multiplication *(moved here from "Analytic Methods")*
- ✅ The Frobenius endomorphism over 𝔽ₚ — trace `a`, `#E(𝔽ₚ) = p+1−a`, the
  Hasse bound, ordinary vs. supersingular, `π` as an element of an imaginary
  quadratic order. *(applet: brute-force point count + trace on the Hasse line)*
- ✅ Lifting Frobenius / Deuring — discussion **rewritten** (2026 rev): what it
  means to lift Frobenius (an analytic model + algebraic eqn over the algebraic
  integers that (1) describes the model, (2) reduces to E/𝔽ₚ mod a prime ideal,
  (3) has an endomorphism descending to φ); Deuring guarantees a lift exists.
  *(applet kept: mult-by-α on Λ = ℤ + αℤ.)*
  - ✅ Worked explicit example `y² = x³ + 3x` over 𝔽₅: #E=10 ⟹ a=−4 ⟹ φ acts as
    α=−2+i; Λ=ℤ[i] forced (order-4 automorphism i); model y²=x³+3x (j=1728)
    descends mod 5; reduce [−2+i] mod 𝔭=(−2+i) (i≡2) ⟹ (x⁵,y⁵). The giant
    lift map x-coordinate shown (\scriptsize, fits) in a collapsed expander.

### §3 — Pictures from lifts of Frobenius  *(NEW — Story 1.5)*
- ✅ "Make a picture using a lift of Frobenius" — the general recipe (ambient
  X(ℂ); plot Fix(F)=E(𝔽ₚ), Fix(Fⁿ)=E(𝔽_{pⁿ}); F = Galois action; orbit length =
  field-of-definition degree) + the **multiplicative group** applet:
  (pⁿ−1)-th roots of unity on the unit circle, coloured by Frobenius-orbit
  length, with a selectable point lighting up its orbit under x↦xᵖ. Point count
  capped (n-range depends on p; p-/n-keyed sliders avoid stale state).
- ⬜ The same idea applied to **elliptic curves** (Fix(αⁿ) on ℂ/Λ = E(𝔽_{pⁿ})).

## Equivalence of Categories

### §4 — Isogenies  *(Story 2)*  *(was §3)*
- ✅ Isogenies: kernels, degree — algebraic vs. analytic descriptions.
  *(concept tab + static Λ ⊆ Λ′ index-3 illustration)*
- ✅ Analytic pictures of isogenies — folding the fundamental domain.
  *(applet: ℓ slabs of ℂ/Λ fold ℓ-to-1 onto ℂ/Λ′; shows the fibre of P)*
- ✅ Vélu's formulas — the algebraic counterpart (demo).
  *(recipe + closed-form 2-isogeny; 𝔽ₚ applet computes codomain E/C and point
  images, verified images land on E/C and #E = #E/C)*

### §5 — Modular tools  *(Story 3)*
- ⬜ The j-function & modular forms.
- ⬜ Modular curves X₀(ℓ).
- ⬜ Modular polynomials Φ_ℓ — reading off isogenies without computing them.

### §6 — The CM bijection  *(Story 4 — closes the loop)*
- ⬜ Quadratic forms, lattice classes & the class-group action.
- ⬜ The bijection: lattice classes ↔ curves in an isogeny class.
- ⬜ Labelling → the bijection.

### Miscellaneous topics  *(catch-all, unnumbered, last)*
- ✅ Isogenies over 𝔽ₚ & the volcano structure — **moved here** from Isogenies
  (discussed too early in the main arc; kept, not deleted).
  *(schematic applet: crater rim + descending ℓ-ary trees, params ℓ/rim/depth)*
- ⬜ Quadratic twists, Hilbert class polynomials, the analytic class number
  formula, supersingular specifics, … — slotted in as they come up.

## Notes
- Two existing tabs move out of the "Analytic Methods" expander: ℂ-tori → §1,
  Endomorphisms/CM → §2. After the move, "Algebraic Curves" / "Analytic Methods"
  expander names are retired in favour of the §1–§6 sections above.
