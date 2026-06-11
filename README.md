# ecfplat

Code supporting computations related to elliptic curves over finite fields via an explicit bijection between lattice classes and elliptic curves in a given isogeny class.

## Overview

Given a pair `(a, p)` with `p` prime and `a² < 4p`, the code works with the isogeny class of elliptic curves over **F**_p whose Frobenius has characteristic polynomial `x² - ax + p`. The central object is an explicit bijection between:

- **Lattice classes** with CM by a root of `x² - ax + p` (equivalently, classes of positive definite binary quadratic forms of discriminant `a² - 4p`)
- **Elliptic curves** in the corresponding isogeny class over **F**_p

This bijection is used as the starting point for computing properties of the elliptic curves (e.g. Mordell–Weil group structure) by working on the lattice side.

## Repository structure

```
notebooks/        # Jupyter notebooks (published)
  userguide.ipynb     # Worked examples and basic use cases

experiments/      # Local scratch notebooks (not tracked by git)

pages/            # Streamlit multi-page app pages
  1_Isogeny_Class.py  # Isogeny class browser (table + lattice picture)
  2_EC_Search.py      # Single-curve lookup (classical + lattice picture)

pycode/           # Core Python library
  alg_classes.py      # General algebraic structures: AbGrp, Ring, RingElement,
                      #   MatrixElement, ZnProduct, Mat_n_Z, Polynomial, PolyFp
  nt.py               # Number theory utilities (gcd, primality, quadratic symbols, …)
  identities.py       # Algebraic identities used in isogeny computations
  qfs.py              # Quadratic form / lattice utilities and modular group action
  modularpolynomials.py  # Atkin and Hilbert modular polynomial data and evaluation
  ecfp.py             # Elliptic curve over F_p computations (j-invariants, models,
                      #   isogeny graphs, Frobenius)
  ecqf_tools.py       # Bijection utilities, Frobenius matrices, Mordell–Weil
                      #   computations, ECQFIsogenyClass
  ecqf.py             # Supporting bijection algorithms (ordinary and supersingular)
  data/               # Precomputed bijection data (JSON)
```

## Web app

A Streamlit app provides a point-and-click interface:

```bash
pip install -r requirements.txt
streamlit run app.py
```

The app has two entry points, accessible from the home page:

- **Isogeny Class** — enter a pair `(a, p)`, load the isogeny class, browse the full bijection table, click a row to see its lattice picture, and download the computed points as a text file.
- **EC Search** — enter coefficients `(f, g, p)` for a curve `y² = x³ + fx + g (mod p)`, look up its trace of Frobenius and associated lattice data, and view both the classical affine point picture and the lattice picture.

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
```

## Precomputed data

The `data/` directory contains precomputed bijections for ordinary pairs `(a, p)` and supersingular primes `p` with `4 ≤ p ≤ 1024`. Use `get_aps_pc()` and `get_ssps_pc()` to list all available pairs.
