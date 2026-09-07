# The removed X₀(ℓ) toolkit — design notes for its successor

The 2026-08 cleanup removed a cluster of Γ₀(ℓ)-side functions from `qfs.py`
(and one baked data table from `modular_viz.py`) that had no remaining
callers.  The full source is one command away (`git show pre-cleanup:pycode/qfs.py`);
this note records *what the code did and how*, so a future X₀(ℓ) module can be
designed cleanly rather than salvaged.

## What existed

All of it operated on binary quadratic forms `(a, b, c)` as proxies for points
τ of the upper half-plane / CM lattices `⟨1, τ⟩`.

- **`qf_to_fun_dom(qf) -> (qf', M)`** — SL₂(ℤ) reduction that *tracked the
  transformation matrix*: same S/T reduction steps as the live `_qf_reduce`,
  but returning the `M ∈ SL₂(ℤ)` with `qf' = M · qf`.  This was the primitive
  everything else built on (you need the matrix, not just the reduced form, to
  talk about Γ₀(ℓ)-cosets).  It went through the generic `MatrixElement`
  classes, which validated every product — the reason the fast fork exists.

- **`find_rrep_g0(M, ℓ)`** — given `M ∈ SL₂(ℤ)`, find the unique right-coset
  representative `M₀` of Γ₀(ℓ)\SL₂(ℤ) (from `gamma_0_coset_reps`, which is
  still live) with `M·M₀⁻¹ ∈ Γ₀(ℓ)` (lower-left entry ≡ 0 mod ℓ).

- **`qf_to_gamma_0_fd(qf, ℓ)` / `qf_mod_gamma_0(qf, ℓ)`** — the composite:
  reduce a form into a chosen fundamental domain for **Γ₀(ℓ)** (SL₂(ℤ)
  fundamental domain translated by the ℓ+1 coset reps), i.e. a canonical
  representative of a point of X₀(ℓ).  This is the Γ₀(ℓ) analog of the live
  `qf_mod_gamma`.

- **`x0_endos_all(p)`** — for every trace `a` with `a² < 4p`, enumerate for
  each form of disc `a²−4p` its Γ₀(p)-orbit points fixed under Fricke
  (via the still-live `qf_x0_endos`), i.e. the points of X₀(p) representing
  endomorphisms — the Eichler/Brandt count behind the 2p+2 identity.

- **`iso_taus_x0_l(qf, ℓ)` / `isos_x0_l_all(d, ℓ)`** — points of X₀(ℓ) lying
  over a given class (resp. all classes of disc d) whose Fricke image lands in
  the same class set: the ℓ-isogeny correspondences between CM points, as
  actual points on the modular curve rather than as graph edges.

- **`modular_viz._X0_ENDOS`** — four small hand-baked tables of the above
  endomorphism points for the applet ℓ's, kept "for the later algebraic
  before/after comparison" that never got built.

## What survives (still live)

`gamma_0_coset_reps(ℓ)`, `gamma_0_orb(qf, ℓ)`, `fricke_inv(qf, ℓ)`, and
`qf_x0_endos(qf, ℓ)` — used by the Background §8 applets.

## Suggestions for the rewrite

- Make the matrix-tracking reduction a *mode* of one reduction function (or
  return the word in S,T and derive the matrix), instead of maintaining two
  parallel implementations again.
- Decide the Γ₀(ℓ) fundamental-domain convention once, document it, and test
  it against `gamma_0_coset_reps` (uniqueness of the representative was
  asserted at runtime in `find_rrep_g0` — make that a unit test instead).
- The consumers to design for: the §8 X₀(ℓ) applets (walker, endomorphism
  picture) and the p–ℓ duality experiments (Eichler counts, Fricke quotients).
