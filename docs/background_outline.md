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
- ✅ The same idea applied to **elliptic curves**: recap (X=ℂ/Λ, F=×α,
  Fix(αⁿ)=E(𝔽_{pⁿ}), count |αⁿ−1|²) + a side-by-side **classical vs lattice**
  applet over the four j=1728 twists y²=x³+ax /𝔽₅ (α an associate of −2+i:
  a=1→1+2i, a=2→2−i, a=3→−2+i, a=4→−1−2i). Lattice picture shows 0 at the
  corner; an 𝔽₂₅ toggle adds Fix(α²) points (the advantage the classical grid
  can't match). **§3 complete.**

## Equivalence of Categories

### §4 — Isogenies  *(Story 2)*  *(was §3)*
- ✅ Isogenies: kernels, degree — algebraic vs. analytic descriptions.
  *(concept tab + static Λ ⊆ Λ′ index-3 illustration)*
- ✅ Analytic pictures of isogenies — folding the fundamental domain.
  *(applet: ℓ slabs of ℂ/Λ fold ℓ-to-1 onto ℂ/Λ′; shows the fibre of P)*
- ✅ Vélu's formulas — the algebraic counterpart (demo).
  *(recipe + closed-form 2-isogeny; 𝔽ₚ applet computes codomain E/C and point
  images, verified images land on E/C and #E = #E/C)*

### §5 — Modular tools  *(Story 3)*  *(built 2026-07, from the MAA talk)*
- ✅ Tab "The j-function and modular curves" — j classifies curves up to iso;
  X(1) = SL₂(ℤ)\ℋ ≅ ℂ; X₀(ℓ) parametrizes cyclic ℓ-isogenies with source/target
  maps to X(1); genus-0 Atkin primes.
- ✅ Tab "Modular polynomials Φ_ℓ" — Φ_ℓ(j(E), j(E′)) = 0 ⟺ cyclic ℓ-isogeny;
  neighbours = roots of Φ_ℓ(j, Y) mod p (no kernels); size/growth ⟹ the
  bootstrapping loop behind 97.04% → 99.86%.

### §6 — The CM bijection  *(Story 4 — closes the loop)*  *(built 2026-07)*
Rebuilt directly on the MAA talk's Act II. Three tabs:
- ✅ "The gallery problem and the equivalence" — gallery problem + Setup (the two
  categories 𝒜(α), ℰ_{𝔽ₚ}(a)) + lifting = equivalence 𝓕_𝔓 + cross-characteristic
  aside.
- ✅ "Isogeny graphs and rigidity" — equivalence ⟹ ℓ-graph isos; **ported disc
  −368 volcano widget** (pills relabel the one graph as lattices or j mod p for
  many (a,p)); rigid spanning set S.
- ✅ "The algorithm and the gallery" — initial values (𝓕(𝒪), orientation) +
  forced propagation (Cayley graphs of Cl(𝒪)); **ported gallery widget**
  (rotate/reflect/reset over the 12 labelings, with mini-torus curve cards).

Graphics ported from the talk via `pycode/slide_viz.py` + `st.components.v1.html`:
the interactive **CM torus** (in §3, "Elliptic curves" tab), the **volcano**, and
the **gallery**. Self-contained canvas/JS with their own validated data.

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
