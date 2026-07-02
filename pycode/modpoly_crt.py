"""Classical modular polynomial Phi_p(X, Y) via the nonstandard CRT method.

Steps 3-4 of the CRT modular-polynomial project, on top of the Hilbert library
(hilbert_crt).  The pipeline:

1. GENERAL FORM.  Phi_p is symmetric, monic of degree p+1 in X, and Phi_p(X,X)
   has degree 2p.  So in the basis X^i Y^j + X^j Y^i the only monomial with
   i = p+1 is X^{p+1} + Y^{p+1} (coefficient 1, known), and the unknowns are
   a_ij for 0 <= j <= i <= p -- (p+1)(p+2)/2 of them.

2. DIAGONAL.  Phi_p(X,X) = - prod H_d(X)^{m_d} over the support: every disc d
   of a complex order containing an element of norm p, i.e. the union of the
   discriminant closures of a^2 - 4p for 0 <= a < 2 sqrt(p); m_d = 1 when
   p | d (d = -4p, or -p for p = 3 mod 4), else 2.  The degree identity
   sum m_d h(d) = 2p is asserted.  This pins the sum of the coefficients on
   each equal-total-degree anti-diagonal, in characteristic 0.

3. INTERPOLATION MOD ell.  For each prime ell != p with precomputed bijections,
   the diagonal relations plus evaluations Phi_p(j1, j2) = 0 at p-isogenous
   pairs over F_ell (read off the CM-lattice side: within a class (a, ell),
   p-neighbour forms = p-isogenous curves) give a linear system for the a_ij
   over F_ell.  Solved by Gaussian elimination; ell is skipped when the system
   is underdetermined, and an inconsistent system raises (it means bad data).

4. CRT.  Coefficients are reconstructed from enough primes; "enough" uses the
   Broker-Sutherland height bound  ln height(Phi_p) <= 6p ln p + 16p + 14 sqrt(p) ln p,
   which is provable but loose (~1.6x on the log scale), so certification is
   conservative.  Output is the dense (p+2)x(p+2) matrix M[i][j] used by
   classical_modpoly / modpoly_nbrs, validated against the Kronecker congruence
   Phi_p = (X^p - Y)(X - Y^p) mod p before returning.
"""

import math
import json
from pathlib import Path

import numpy as np

from nt import primeQ, crt_list
from identities import disc_closure
from qfs import get_qfs_strict, qf_isogs_hor
from ecqf_tools import ecqf_ord_1K_pc, ecqf_ss_1K_pc
from alg_classes import poly_ring, ZZ
from modularpolynomials import hilb_polys_dict

_DATA_DIR = Path(__file__).parent / 'data'


from hilbert_crt import hilbert_library


########################
# Step 2: the diagonal #
########################

def diagonal_support(p: int) -> dict:
    """{d: multiplicity} with Phi_p(X,X) = - prod H_d^{m_d}: the discs of complex
    orders containing a norm-p element; m_d = 1 when p | d, else 2."""
    supp = {}
    a = 0
    while a * a < 4 * p:
        for d in disc_closure(a * a - 4 * p):
            supp[d] = 1 if d % p == 0 else 2
        a += 1
    deg = sum(m * len(get_qfs_strict(d)) for d, m in supp.items())
    if deg != 2 * p:
        raise ValueError(f'diagonal degree check failed for p={p}: {deg} != {2*p}')
    return supp


def phi_diagonal(p: int, hilb: dict = None):
    """Phi_p(X, X) in characteristic 0, as a Poly over ZZ (leading coefficient -1)."""
    if hilb is None:
        hilb = hilbert_library()
    supp = diagonal_support(p)
    missing = sorted(d for d in supp if d not in hilb)
    if missing:
        raise ValueError(f'Hilbert polynomials missing for {missing}; '
                         f'run hilbert_crt.hilbert_poly_search on them first')
    Zx = poly_ring(ZZ)
    F = Zx.one
    for d, m in sorted(supp.items()):
        F = F * Zx.poly(hilb[d]) ** m
    return -1 * F


##################################
# Step 3: solve Phi_p mod ell    #
##################################

def phi_monomials(p: int) -> list[tuple]:
    """The unknown symmetric monomials (i, j), 0 <= j <= i <= p."""
    return [(i, j) for i in range(p + 1) for j in range(i + 1)]


def isogenous_pairs_mod_ell(p: int, ell: int) -> list[tuple]:
    """Pairs (j1, j2) with Phi_p(j1, j2) = 0 over F_ell, from the precomputed
    ordinary bijections: within a class (a, ell), the p-neighbour forms of a
    curve's form are its p-isogenous curves (both orientations give valid pairs,
    so the global conjugation ambiguity is harmless)."""
    pairs = set()
    for (a, l0), bij in ecqf_ord_1K_pc.items():
        if l0 != ell or a <= 0:
            continue
        qf_to_j = {qf: j for j, qf in bij.items()}
        for j1, qf in bij.items():
            for qf2 in qf_isogs_hor(qf, p):
                j2 = qf_to_j.get(tuple(qf2))
                if j2 is not None:
                    pairs.add((min(j1, j2), max(j1, j2)))
    ss = ecqf_ss_1K_pc.get(ell, {})            # supersingular class: oriented CM under
    qf_to_sig = {qf: sig for sig, qf in ss.items()}       # an order of Q(sqrt(-ell))
    for (j1, s1), qf in ss.items():
        for qf2 in qf_isogs_hor(qf, p):
            sig2 = qf_to_sig.get(tuple(qf2))
            if sig2 is not None:
                pairs.add((min(j1, sig2[0]), max(j1, sig2[0])))
    return sorted(pairs)


def special_value_rows(p: int, hilb: dict = None) -> dict:
    """Char-0 linear relations on the a_ij from the class-number-1 j0's.

    For h(d0) = 1 the p+1 cyclic p-isogenies from E_{j0} either stay at disc d0
    (their codomain is E_{j0} itself, 1 + chi of them, chi = kronecker(d0, p))
    or descend to disc d0 p^2, collapsing under Aut(E_{j0}) into the h(d0 p^2)
    roots of H_{d0 p^2} with multiplicity u = 3, 2, 1 for d0 = -3, -4, other.
    In the standard normalization everything is monic:

        Phi_p(X, j0) = (X - j0)^{1+chi} * H_{d0 p^2}(X)^u.

    Each X^i coefficient (i <= p) is one linear relation, valid modulo every ell.
    Only the j0 with H_{d0 p^2} in the library contribute.  Returns
    {'rows': dense rows over phi_monomials(p), 'rhs': ints, 'used': [d0...]}."""
    from modularpolynomials import heeg_js
    from nt import quad_rec
    if hilb is None:
        hilb = hilbert_library()
    mons = phi_monomials(p)
    idx = {m: k for k, m in enumerate(mons)}
    rows, rhs, used = [], [], []
    Zx = poly_ring(ZZ)
    for d0 in sorted(heeg_js, reverse=True):
        dd = d0 * p * p
        if dd not in hilb:
            continue
        j0 = heeg_js[d0]
        chi = quad_rec(d0 % p, p)
        u = {-3: 3, -4: 2}.get(d0, 1)
        F = (Zx.poly([-j0, 1]) ** (1 + chi)) * (Zx.poly(hilb[dd]) ** u)
        cs = F.int_coefs() + [0] * (p + 2)
        j0pow = [j0 ** k for k in range(p + 2)]
        for i in range(p + 1):                 # i = p+1 is the trivial identity 1 = 1
            row = [0] * len(mons)
            for j in range(p + 1):
                row[idx[(max(i, j), min(i, j))]] += j0pow[j]
            rows.append(row)
            rhs.append(cs[i] - (j0pow[p + 1] if i == 0 else 0))
        used.append(d0)
    return {'rows': rows, 'rhs': rhs, 'used': used}


def _solve_unique_mod(rows: list, rhs: list, n: int, ell: int):
    """The unique solution of the linear system over F_ell, or None when the rank
    is < n.  Raises on an inconsistent system -- that means bad isogeny data."""
    A = np.array([[c % ell for c in r] + [b % ell] for r, b in zip(rows, rhs)],
                 dtype=np.int64)               # reduce in Python first: rhs holds char-0 ints
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, len(A)) if A[i, c] % ell), None)
        if piv is None:
            return None                        # rank-deficient in column c
        A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), -1, ell)) % ell
        mask = np.arange(len(A)) != r
        A[mask] = (A[mask] - np.outer(A[mask, c], A[r])) % ell
        r += 1
    if np.any(A[r:, n] % ell):
        raise ValueError(f'inconsistent linear system mod {ell} -- bad isogeny pair data')
    return [int(A[i, n]) for i in range(n)]


def solve_phi_from_pairs(p: int, ell: int, diag_coefs: list, pairs: list,
                         extra: dict = None) -> dict:
    """{(i, j): a_ij mod ell} for the unknown coefficients of Phi_p, or None when
    the equations available at ell leave the system underdetermined.

    Equations: (a) for each k <= 2p, the anti-diagonal sum equals the char-0
    diagonal coefficient of X^k (minus the known X^{p+1}+Y^{p+1} contribution at
    k = p+1); (b) for each supplied p-isogenous pair, Phi_p(j1, j2) = 0;
    (c) optional char-0 rows from special_value_rows, shared across all ell."""
    mons = phi_monomials(p)
    idx = {m: k for k, m in enumerate(mons)}
    rows, rhs = [], []
    if extra is not None:
        rows += extra['rows']
        rhs += extra['rhs']
    for k in range(2 * p + 1):                 # (a) anti-diagonal sums
        row = [0] * len(mons)
        for (i, j), t in idx.items():
            if i + j == k:
                row[t] = 1 if i == j else 2
        rows.append(row)
        rhs.append((diag_coefs[k] if k < len(diag_coefs) else 0) - (2 if k == p + 1 else 0))
    for j1, j2 in pairs:                       # (b) evaluations
        pw1 = [pow(j1, i, ell) for i in range(p + 2)]
        pw2 = [pow(j2, i, ell) for i in range(p + 2)]
        row = [0] * len(mons)
        for (i, j), t in idx.items():
            row[t] = (pw1[i] * pw2[j] + (pw1[j] * pw2[i] if i != j else 0)) % ell
        rows.append(row)
        rhs.append(-(pw1[p + 1] + pw2[p + 1]))
    sol = _solve_unique_mod(rows, rhs, len(mons), ell)
    if sol is None:
        return None
    return {m: sol[t] for m, t in idx.items()}


def solve_phi_mod_ell(p: int, ell: int, diag_coefs: list, extra: dict = None) -> dict:
    """solve_phi_from_pairs with the pairs read off the precomputed (ell < 1024)
    bijections."""
    return solve_phi_from_pairs(p, ell, diag_coefs, isogenous_pairs_mod_ell(p, ell), extra)


##############################
# Step 4: CRT + certification #
##############################

def modpoly_height_bound_log2(p: int) -> float:
    """Broker-Sutherland: ln height(Phi_p) <= 6p ln p + 16p + 14 sqrt(p) ln p."""
    return (6 * p * math.log(p) + 16 * p + 14 * math.sqrt(p) * math.log(p)) / math.log(2)


def phi_p_via_crt(p: int, ells: list = None, verbose: bool = True) -> dict:
    """Phi_p(X, Y) by CRT interpolation across the precomputed primes.

    Returns {'M': dense (p+2)x(p+2) matrix, 'primes': ells used, 'certified': ...};
    M is in the classical_modpoly format and is checked against the Kronecker
    congruence mod p before returning."""
    if not primeQ(p):
        raise ValueError(f'p = {p} must be prime')
    hilb = hilbert_library()
    diag = phi_diagonal(p, hilb).int_coefs()
    extra = special_value_rows(p, hilb)
    if verbose and extra['used']:
        print(f'  special-value relations from j0 of disc {extra["used"]}')
    if ells is None:
        ells = sorted({l0 for a, l0 in ecqf_ord_1K_pc if a > 0}, reverse=True)
    bound = modpoly_height_bound_log2(p)
    residues, bits, skipped = [], 0.0, []
    for ell in ells:
        if ell == p:
            continue
        if bits > bound + 1:
            break
        sol = solve_phi_mod_ell(p, ell, diag, extra)
        if sol is None:
            skipped.append(ell)
            continue
        residues.append((sol, ell))
        bits += math.log2(ell)
        if verbose and len(residues) % 10 == 0:
            print(f'  {len(residues)} primes, {bits:.0f} / {bound + 1:.0f} bits', flush=True)
    if not residues:
        return {'p': p, 'status': 'fail', 'reason': 'no prime gave a determined system',
                'skipped': skipped}
    M = [[0] * (p + 2) for _ in range(p + 2)]
    M[p + 1][0] = M[0][p + 1] = 1
    for m in phi_monomials(p):
        c, mod = crt_list([(sol[m], ell) for sol, ell in residues])
        M[m[0]][m[1]] = M[m[1]][m[0]] = c
    rec = {'p': p, 'status': 'ok', 'M': M, 'primes': [ell for _, ell in residues],
           'log2_modulus': bits, 'log2_bound': bound, 'certified': bits > bound + 1,
           'skipped': skipped, 'kronecker_ok': kronecker_check(p, M)}
    if not rec['kronecker_ok']:
        rec['status'] = 'suspect'
    return rec


##############
# Validation #
##############

def kronecker_check(p: int, M: list) -> bool:
    """Phi_p(X, Y) = (X^p - Y)(X - Y^p) mod p, i.e. X^{p+1} + Y^{p+1} - X^p Y^p - XY."""
    want = {(p + 1, 0): 1, (0, p + 1): 1, (p, p): -1, (1, 1): -1}
    for i in range(p + 2):
        for j in range(p + 2):
            if (M[i][j] - want.get((i, j), 0)) % p:
                return False
    return True


def check_phi_against_cache(p: int, M: list) -> bool:
    """Exact comparison with the validated q-expansion computation (l <= 23 cached)."""
    from modularpolynomials import classical_modpoly
    return M == classical_modpoly(p)
