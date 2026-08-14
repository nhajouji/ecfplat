"""Derivation pipeline for algebraic models of modular curves (paper appendix A).

Step 0 of the model-library project: the universal Tate curve

    E_{u,v}:  y^2 + (1-u)*x*y - v*y = x^3 - v*x^2,   P0 = (0,0)

over Spec Q[u,v], multiples [m]P0 as rational functions in (u,v), and the
planar curves C_ell = {r_ell(u,v) = 0} that are birational to X_1(ell).

Everything here is exact (sympy over QQ).  Run as a script for the
verification battery: r_5 = u - v, r_7 = u^3 + u*v - v^2 (paper values),
plus r_11 vanishing on the appendix model of X_1(11).
"""
import sympy as sp

u, v = sp.symbols('u v')

# Tate normal form coefficients: y^2 + a1 xy + a3 y = x^3 + a2 x^2 (a4=a6=0)
A1, A2, A3 = 1 - u, -v, -v


def _neg(P):
    x, y = P
    return (x, -y - A1 * x - A3)


def _add(P, Q):
    """Group law on E_{u,v}; P, Q distinct non-opposite affine points."""
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and y1 == y2:
        lam = sp.cancel((3 * x1**2 + 2 * A2 * x1 - A1 * y1) / (2 * y1 + A1 * x1 + A3))
    else:
        lam = sp.cancel((y2 - y1) / (x2 - x1))
    nu = sp.cancel(y1 - lam * x1)
    x3 = sp.cancel(lam**2 + A1 * lam - A2 - x1 - x2)
    y3 = sp.cancel(-(lam * x3 + nu) - A1 * x3 - A3)
    return (x3, y3)


def multiples(n):
    """[1]P0 .. [n]P0 as pairs of rational functions in (u, v)."""
    P0 = (sp.Integer(0), sp.Integer(0))
    out = [P0]
    P = P0
    for _ in range(n - 1):
        P = _add(P, P0) if P != P0 else _add(P0, P0)
        out.append(P)
    return out


def j_universal():
    """j-invariant of E_{u,v} as a rational function in (u, v)."""
    b2 = A1**2 + 4 * A2
    b4 = A1 * A3
    b6 = A3**2
    b8 = A2 * A3**2  # a4 = a6 = 0 case of the general b8
    c4 = b2**2 - 24 * b4
    disc = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    return sp.cancel(c4**3 / disc)


def r_ell(ell, mults=None):
    """The polynomial r_ell(u, v) cutting out C_ell (P0 has exact order ell).

    Sets x([a]P0) = x([b]P0) with a+b = ell adjacent, then strips the factors
    coming from lower-order degenerations: a factor F is spurious iff it also
    divides some x([i]P0) = x([k]P0) condition with i+k < ell (or a
    denominator, i.e. a smaller [n]P0 = O locus).  We identify r_ell as the
    product of factors that survive.
    """
    a, b = (ell - 1) // 2, (ell + 1) // 2
    M = mults if mults is not None else multiples(b)
    diff = sp.together(M[a - 1][0] - M[b - 1][0])
    num = sp.factor(sp.numer(diff))
    # denominators of all multiples = smaller-torsion loci (spurious)
    spurious = sp.Integer(1)
    for (x, y) in M[:b - 1]:
        spurious *= sp.denom(sp.together(x))
    # x-coordinate collisions of strictly smaller total order
    for i in range(1, b):
        for k in range(i, b):
            if i + k < ell and (i, k) != (a, b):
                d = sp.together(M[i - 1][0] - M[k - 1][0]) if i != k else sp.Integer(1)
                if d != 1:
                    spurious *= sp.numer(d)
    keep = sp.Integer(1)
    for fac, _mult in sp.factor_list(num)[1]:
        if not fac.is_number and not _divides(fac, spurious):
            keep *= fac
    return sp.expand(keep)


def _divides(f, g):
    """True if polynomial f divides polynomial g in Q[u,v]."""
    q = sp.cancel(g / f)
    return sp.denom(sp.together(q)) == 1


# ---------------------------------------------------------------------------
# Fast path: division polynomials evaluated at P0 = (0,0).
#
# psi_n(P0) in Z[u,v] via the standard doubling recurrences
#   psi_{2n+1} = psi_{n+2} psi_n^3 - psi_{n-1} psi_{n+1}^3
#   psi_{2n}   = psi_n (psi_{n+2} psi_{n-1}^2 - psi_{n-2} psi_{n+1}^2) / psi_2
# which need psi at a *general* point... but evaluated along the multiples of
# P0 they close up if we use the elliptic-net (EDS) form: the sequence
# W(n) = psi_n(P0) is a proper divisibility sequence satisfying
#   W(m+n) W(m-n) W(r)^2 = W(m+r) W(m-r) W(n)^2 - W(n+r) W(n-r) W(m)^2,
# and the doubling formulas above hold verbatim for W.  Initial values
# W(1)..W(4) are computed from the curve coefficients at P0.
# ---------------------------------------------------------------------------

def _eds_initial():
    """W(1)..W(4) = psi_1..psi_4 at P0 = (0,0), as Poly in ZZ[u,v]."""
    x, y = sp.symbols('x y')
    b2 = A1**2 + 4 * A2
    b4 = A1 * A3
    b6 = A3**2
    b8 = A2 * A3**2
    psi2 = 2 * y + A1 * x + A3
    psi3 = 3 * x**4 + b2 * x**3 + 3 * b4 * x**2 + 3 * b6 * x + b8
    psi4 = psi2 * (2 * x**6 + b2 * x**5 + 5 * b4 * x**4 + 10 * b6 * x**3
                   + 10 * b8 * x**2 + (b2 * b8 - b4 * b6) * x + (b4 * b8 - b6**2))
    at0 = {x: 0, y: 0}
    return [sp.Poly(sp.expand(e.subs(at0)), u, v) for e in
            (sp.Integer(1), psi2, psi3, psi4)]


def eds_sequence(nmax):
    """W(1)..W(nmax): division polynomial values at P0, exact in ZZ[u,v]."""
    W = {1: None, 2: None, 3: None, 4: None}
    W[1], W[2], W[3], W[4] = _eds_initial()
    one = sp.Poly(1, u, v)
    W[0] = sp.Poly(0, u, v)

    def get(n):
        if n in W:
            return W[n]
        if n % 2:
            m = (n - 1) // 2
            val = get(m + 2) * get(m)**3 - get(m - 1) * get(m + 1)**3
        else:
            m = n // 2
            val = (get(m) * (get(m + 2) * get(m - 1)**2
                             - get(m - 2) * get(m + 1)**2)).exquo(W[2])
        W[n] = val
        return val

    return [get(n) for n in range(1, nmax + 1)]


def r_ell_fast(ell, Wseq=None):
    """r_ell(u,v) from W(ell) = psi_ell(P0): strip factors dividing any
    smaller W(n) (lower-order torsion / degenerate loci)."""
    W = Wseq if Wseq is not None else eds_sequence(ell)
    target = W[ell - 1]
    keep = sp.Poly(1, u, v)
    for fac, mult in sp.factor_list(target.as_expr())[1]:
        f = sp.Poly(fac, u, v)
        if any(Wn.rem(f).is_zero for Wn in W[:ell - 1] if not Wn.is_zero):
            continue
        keep = keep * f**mult
    return keep.as_expr()
