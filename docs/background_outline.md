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

### §1 — Elliptic curves: the basics  *(foundation / review)*
- ✅ Algebraic curves over ℝ
- ✅ Real elliptic curves + group law
- ✅ Elliptic curves over 𝔽ₚ + group law
- ✅↪ Elliptic curves over ℂ: complex tori *(moved here from "Analytic Methods")*

### §2 — Frobenius and the lattice pictures  *(Story 1)*
- ✅↪ Endomorphisms & Complex Multiplication *(moved here from "Analytic Methods")*
- ✅ The Frobenius endomorphism over 𝔽ₚ — trace `a`, `#E(𝔽ₚ) = p+1−a`, the
  Hasse bound, ordinary vs. supersingular, `π` as an element of an imaginary
  quadratic order. *(applet: brute-force point count + trace on the Hasse line)*
- ✅ Lifting Frobenius / Deuring — the Deuring lifting theorem; the lattice plus
  the marked lift of `π` *is* the picture we display. Closes Story 1.
  *(applet: mult-by-π on Λ = ℤ + πℤ, fundamental cell → rotated/scaled image)*

### §3 — Isogenies  *(Story 2)*
- ✅ Isogenies: kernels, degree — algebraic vs. analytic descriptions.
  *(concept tab + static Λ ⊆ Λ′ index-3 illustration)*
- ✅ Analytic pictures of isogenies — folding the fundamental domain.
  *(applet: ℓ slabs of ℂ/Λ fold ℓ-to-1 onto ℂ/Λ′; shows the fibre of P)*
- ✅ Vélu's formulas — the algebraic counterpart (demo).
  *(recipe + closed-form 2-isogeny; 𝔽ₚ applet computes codomain E/C and point
  images, verified images land on E/C and #E = #E/C)*
- ✅ Isogenies over 𝔽ₚ & the volcano structure.
  *(schematic applet: crater rim + descending ℓ-ary trees, params ℓ/rim/depth)*

### §4 — Modular tools  *(Story 3)*
- ⬜ The j-function & modular forms.
- ⬜ Modular curves X₀(ℓ).
- ⬜ Modular polynomials Φ_ℓ — reading off isogenies without computing them.

### §5 — The CM bijection  *(Story 4 — closes the loop)*
- ⬜ Quadratic forms, lattice classes & the class-group action.
- ⬜ The bijection: lattice classes ↔ curves in an isogeny class.
- ⬜ Labelling → the bijection.

### §6 — Miscellaneous topics  *(catch-all)*
- ⬜ Quadratic twists, Hilbert class polynomials, the analytic class number
  formula, supersingular specifics, … — slotted in as they come up.

## Notes
- Two existing tabs move out of the "Analytic Methods" expander: ℂ-tori → §1,
  Endomorphisms/CM → §2. After the move, "Algebraic Curves" / "Analytic Methods"
  expander names are retired in favour of the §1–§6 sections above.
