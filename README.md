# ecfplat

Code supporting computations related to elliptic curves over finite fields via an explicit bijection between lattice classes and elliptic curves in a given isogeny class.

The lattice point data produced here is used as input for shader-rendered artwork by Nadir Hajouji and Steve Trettel, displayed at [elliptic-curves.art](https://elliptic-curves.art/).

## Overview

Given a pair `(a, p)` with `p` prime and `a² < 4p`, the code works with the isogeny class of elliptic curves over **F**_p whose Frobenius has characteristic polynomial `x² - ax + p`. The central object is an explicit bijection between:

- **Lattice classes** with CM by a root of `x² - ax + p` (equivalently, classes of positive definite binary quadratic forms of discriminant `a² - 4p`)
- **Elliptic curves** in the corresponding isogeny class over **F**_p

This bijection is used as the starting point for computing properties of the elliptic curves (e.g. Mordell–Weil group structure) by working on the lattice side.

### How the bijection is computed

On the lattice side the class group acts on quadratic forms by `ℓ`-isogenies. A **rigid l-set** is a small set of primes `ℓ` whose isogeny directions form an *independent generating set* of the class group — together, when there are two or more directions of order > 2, with a "sum" element that pins the relative orientation of the cycles. Walking these directions assigns every class an integer-tuple coordinate `(x₁, …, xₙ)`. The same coordinate labelling is built on the curve side from `F_p` `ℓ`-isogeny adjacency (read off Atkin modular polynomials, no isogeny is computed explicitly), and matching the two labellings coordinate-by-coordinate yields the bijection.

The search for a rigid l-set is `disc_rigid_lset_search` in `pycode/ecqf_bij.py`; the end-to-end ordinary bijection is `ecqf_full_bijection_ord`. The labelling/cube machinery lives in `pycode/graph_tools.py`.

## Repository structure

```
notebooks/        # Jupyter notebooks (published)
  userguide.ipynb     # Worked examples and basic use cases

experiments/      # Local scratch notebooks (not tracked by git)

pages/            # Streamlit multi-page app pages
  0_Homepage.py       # Landing page with project description and navigation
  1_Isogeny_Class.py  # Isogeny class browser: bijection table, isogeny graph,
                      #   and cross-navigation to EC Search
  2_EC_Search.py      # Single-curve lookup: classical and lattice pictures,
                      #   point download, cross-navigation to Isogeny Class
  3_Background.py     # Crash course: interactive lecture tabs covering elliptic
                      #   curves over ℝ/𝔽ₚ/ℂ, isogenies, CM, and Frobenius

pycode/           # Core Python library
  alg_classes.py      # Algebraic structures: AbGrp/Ring/Field + elements, matrix
                      #   rings (Mat_n_Z), prime/extension fields (GF_p, GF_pn),
                      #   Polynomial/PolyFp
  nt.py               # Number theory utilities (gcd, primality, quadratic symbols,
                      #   Frobenius-eigenvalue / isogeny-kernel extension degrees, …)
  identities.py       # Algebraic identities used in isogeny computations
  qfs.py              # Quadratic form / lattice utilities and modular group action
  modularpolynomials.py  # Atkin and Hilbert modular polynomial data and evaluation
  ecfp.py             # Elliptic curves over F_p (j-invariants, models, isogeny
                      #   graphs, Frobenius)
  ecqf_bij.py         # Rigid l-set search and the ordinary lattice <-> curve
                      #   bijection (disc_rigid_lset_search, ecqf_full_bijection_ord)
  rigid_cache.py      # Per-discriminant cache (search + lattice-side bijection) with
                      #   a populate/update CLI and a cached (a, p) entrypoint
  ldata_cache.py      # Per-discriminant rigid-l-set-data cache + populate/update CLI
  ecqf_tools.py       # Bijection utilities, Frobenius matrices, Mordell–Weil
                      #   computations, ECQFIsogenyClass, precomputed-table loaders
  ecqf.py             # Legacy utilities and precomputed-table loading (the bijection
                      #   algorithms now live in ecqf_bij.py)
  graph_tools.py      # Isogeny graph utilities: adjacency matrices, cycle/tree
                      #   decompositions, the Zⁿ labelling algorithm
  graphic_tools.py    # Helpers for the lattice-point artwork output
  misctools.py        # Small shared utilities (dict composition, tuple helpers)
  data/               # Precomputed JSON (see "Precomputed data")
```

## Web app

A Streamlit app provides a point-and-click interface:

```bash
pip install -r requirements.txt
streamlit run app.py
```

The app has four pages:

- **Homepage** — project description and links to the two main tools.
- **Background** — crash course on the underlying mathematics, with interactive applets. Currently implemented: chord-tangent group law on elliptic curves over ℝ, and a τ explorer for the complex lattice ℂ/Λ.
- **EC Search** — enter coefficients `(f, g, p)` for a curve `y² = x³ + fx + g (mod p)`, look up its trace of Frobenius and associated lattice data, view classical and lattice pictures, and navigate directly to its isogeny class.
- **Isogeny Class** — enter a pair `(a, p)`, browse the full bijection table, view degree-ℓ isogeny graphs (adjacency matrix + concentric-ring picture with horizontal/vertical edges distinguished by colour), and navigate to any individual curve in EC Search.

## Getting started (library)

Install dependencies:

```bash
pip install -r requirements.txt
```

Open the user guide:

```bash
jupyter notebook notebooks/userguide.ipynb
```

The notebook walks through:
- Checking whether a pair `(a, p)` has precomputed data
- Initializing an `ECQFIsogenyClass` object
- Viewing the bijection as a pandas DataFrame
- Visualizing lattice classes in the upper half-plane
- Computing Mordell–Weil group data from the lattice side

### Quick example

```python
import sys
sys.path.insert(0, 'pycode/')
from ecqf_tools import ECQFIsogenyClass, ap_in_pc_data

# Check if (a=22, p=1021) is in the precomputed data
ap_in_pc_data((22, 1021))   # True

# Load the isogeny class
isoclass = ECQFIsogenyClass(22, 1021)

# View all lattice classes and their corresponding elliptic curve data
isoclass.ecqf_df()

# Compute the degree-5 isogeny graph adjacency matrix
isoclass.adjacency_matrix(5)
```

## Precomputed data

`pycode/data/` holds two layers of precomputed results, plus supporting modular-polynomial tables.

**Per-`(a, p)` bijections** — the end product, the full curve ↔ lattice-class bijection for each pair:
- `ecqf_ord_pcbij_4_1024.json` — ordinary classes, **6 725** pairs `(a, p)` with `4 ≤ p ≤ 1024`.
- `ecqf_ss_pcbij_4_1024_INC.json` — supersingular classes, keyed by prime `p` (158 keys).
- `ecqf_ord_pcbij_4to256.json` — an earlier / smaller-range variant.

List available keys with `get_aps_pc()` / `get_ssps_pc()` and test membership with `ap_in_pc_data((a, p))` (all in `ecqf_tools.py`).

**Per-discriminant data** — the `(a, p)`-independent layer. The expensive lattice-side computation depends only on `d = a² − 4p`, so it is precomputed once per discriminant for **all 2 048** discriminants in `[−4096, −3]`:
- `qf_ldata.json` — the rigid l-set search result for each `d` (generating primes, their orders, the sum/pinning element). **2 036** of the 2 048 admit a rigid spanning set from the 15 supersingular primes `{2, …, 71}`; **14** of those need a prime-power generator `(ℓ, k)`; the remaining **12** have no rigid l-set in that prime pool.
- `rigid_lset_cache.json` — the same search data, plus the `(a, p)`-independent lattice-side labelling `(x₁, …, xₙ) ↦ (a, b, c)` for each `d`.

Both are regenerated/extended with incremental command-line tools (re-running only fills what is missing):

```bash
python pycode/rigid_cache.py --min -8192   # search + lattice-side bijection
python pycode/ldata_cache.py --min -8192   # search data only (smaller)
```

**Modular-polynomial tables** — used to read `F_p` isogeny adjacency without computing isogenies:
- `atkinpolys.json` — Atkin modular polynomials for the 15 primes `ℓ ∈ {2, …, 71}` whose Atkin–Lehner quotient `X₀(ℓ)⁺` has genus 0 (the supersingular primes).
- `hilbpolys.json`, `jcoefs.json` — Hilbert class polynomials and related data for small discriminants.

### Validation

The per-discriminant data has been checked as follows:

- **Determinism / cache integrity.** The lattice-side labelling is a deterministic function of `d`, so a cached entry reproduces a fresh recomputation *exactly*; verified across the cached range and through a JSON round-trip.
- **Bijectivity.** Every successful search result yields a *complete* bijection of the class group (injective and onto), confirmed by reconstructing the lattice-side labelling for all 2 036 solved discriminants.
- **Regression against the trusted tables.** Rebuilding each `(a, p)` bijection from the per-discriminant cache reproduces the stored `ecqf_ord_pcbij_4_1024.json` for **3 913** of the 6 725 pairs and its exact complex-conjugate for the other **2 812** — **0 disagreements**. (Both are valid bijections; the labelling is canonicalized to a fixed orientation, so the choice is consistent across pairs.)
- **Backward compatibility.** Extending the search to allow prime-power generators recovered the 14 additional discriminants noted above while leaving every prime-only result byte-identical (checked on a 300-discriminant sample).
