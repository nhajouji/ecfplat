## Data loading

import json
from pathlib import Path
from alg_classes import Polynomial
from nt import discfac
from qfs import class_group_id

_DATA_DIR = Path(__file__).parent / 'data'

with open(_DATA_DIR / 'atkinpolys.json', 'r') as f:
    atkin_polys_raw = json.load(f)
atkin_polys_dict = {int(p):atkin_polys_raw[p] for p in atkin_polys_raw}

with open(_DATA_DIR / 'hilbpolys.json', 'r') as f:
    hilbpolys_raw = json.load(f)
hilb_polys_dict = {int(d):hilbpolys_raw[d] for d in hilbpolys_raw}
heeg_js = {d:-hilb_polys_dict[d][0] for d in hilb_polys_dict if len(hilb_polys_dict[d])==2}

heeg_js_to_qfs_dics = {}
for d in heeg_js:
    d0 = discfac(d)[0]
    if d0 == d:
        heeg_js_to_qfs_dics[d]={heeg_js[d]:class_group_id(d)}
    else:
        heeg_js_to_qfs_dics[d] = {heeg_js[d0]:class_group_id(d0),heeg_js[d]:class_group_id(d)}

def small_bij_check(d:int)->dict:
    if d in heeg_js_to_qfs_dics:
        return heeg_js_to_qfs_dics[d]
    else:
        return {}

## Polynomials and Rational Functions

def poly_eval_mod(coefs_lc_to_const:list[int],x:int,p:int,rev = False):
    if rev == True:
        coefs_lc_to_const = coefs_lc_to_const[::-1]
    evx = coefs_lc_to_const[0]%p
    for c in coefs_lc_to_const[1:]:
        evx = (evx*x+c)%p
    return evx


def rat_eval_mod(coefdic_lc_to_const:dict,x:int,p:int,rev = False):
    ev = 1
    for e in coefdic_lc_to_const:
        fac_coefs = coefdic_lc_to_const[e]
        fac_ev = poly_eval_mod(coefs_lc_to_const=fac_coefs,x=x,p=p,rev=rev)
        if e>= 0 or fac_ev %p != 0:
            ev = (ev*pow(fac_ev,e,p)) %p
        else:
            return 'Infinity'
    return ev


## Count roots

def count_roots_fp_bf(coefs_lc_to_const:list[int],p:int)->int:
    return len([x for x in range(p) if poly_eval_mod(coefs_lc_to_const,x,p)==0])

## Evaluating modular polynomials

def eval_atk(x:int,l:int,p:int):
    if l not in atkin_polys_dict:
        raise ValueError('l does not divide order of the Monster group')
    cs = atkin_polys_dict[l]
    a,b1,b3= tuple([poly_eval_mod(ci[::-1],x,p) for ci in cs])
    b = (b1*pow(b3,3,p))%p
    return ((-a)%p,b)


def atk_poly_a(l:int,p=0):
    if l not in atkin_polys_dict:
        raise ValueError('l does not divide order of the Monster group')
    return Polynomial(atkin_polys_dict[l][0],char=p)

#This, in isolation, captures endorphisms of degree l for j = 0
def atk_poly_b1(l:int,p=0):
    if l not in atkin_polys_dict:
        raise ValueError('l does not divide order of the Monster group')
    return Polynomial(atkin_polys_dict[l][1],char=p)

def atk_poly_b3(l:int,p=0):
    if l not in atkin_polys_dict:
        raise ValueError('l does not divide order of the Monster group')
    return Polynomial(atkin_polys_dict[l][2],char=p)

def atk_poly_b(l:int,p=0):
    b1 = atk_poly_b1(l,p)
    b3 = atk_poly_b3(l,p)
    return b1*(b3**3)

def atk_poly_ab(l:int,p=0):
    return (atk_poly_a(l,p),atk_poly_b(l,p))

def atk_at_j(j:int,l:int,p=0):
    j2 = Polynomial([j**2],char=p)
    a,b = atk_poly_ab(l,p)
    return j2-j*a+b


def atk_at_j_fpfac(j:int,l:int,p:int):
    poly = atk_at_j(j,l,p).mod(p)
    return poly.fp_factor()


## Classical modular polynomials Phi_l(X, Y), computed from the j-function q-expansion
# The Atkin polynomials above only cover the 15 Atkin primes (those dividing |Monster|,
# i.e. with X_0(l)^+ of genus 0).  For an arbitrary prime l
# we compute the classical modular polynomial Phi_l(X, Y) directly from the q-expansion of
# the j-function (data/jq_coeffs.json, the coefficients c(n) of j = q^-1 + 744 + ...).
# Phi_l(j, Y) mod p factors over F_p into the l+1 j-invariants l-isogenous to j -- the
# general-l analogue of eval_atk, with no isogeny computed over an extension field.

with open(_DATA_DIR / 'jq_coeffs.json') as f:
    _jq = json.load(f)
_JQ_MIN = _jq['min_n']                                   # -1
_JQ = _jq['coeffs']                                      # _JQ[i] = c(i + _JQ_MIN)
_JQ_NMAX = _JQ_MIN + len(_JQ) - 1
def _jc(n):                                              # j-coefficient c(n), n >= -1
    return _JQ[n - _JQ_MIN]

# integer Laurent series as a dict {exponent: coefficient}
def _smul(A, s): return {e: v*s for e, v in A.items() if v*s}
def _sadd(*Ds):
    R = {}
    for D in Ds:
        for e, v in D.items(): R[e] = R.get(e, 0) + v
    return {e: v for e, v in R.items() if v}
def _smulser(A, B, emax):
    R = {}
    for e1, v1 in A.items():
        for e2, v2 in B.items():
            e = e1 + e2
            if e <= emax: R[e] = R.get(e, 0) + v1*v2
    return {e: v for e, v in R.items() if v}
def _spow(A, k, emax):
    R = {0: 1}
    for _ in range(k): R = _smulser(R, A, emax)
    return R

def compute_modpoly(l:int):
    """Classical modular polynomial Phi_l(X, Y) as a dense (l+2)x(l+2) integer matrix M with
    M[i][j] = coefficient of X^i Y^j, computed from the q-expansion of j.

    Cyclotomic-free method: the l "small" conjugates j((tau+b)/l) have power sums
        P_k = sum_b j((tau+b)/l)^k = l * (l-divisible q^{1/l}-part of j(tau/l)^k)
    (the sum over b kills the non-l-divisible exponents, so no zeta_l is needed), and these
    have shallow poles.  Newton's identities give their elementary symmetric functions
    sigma_k; then Phi_l(X, j) = (X - j(l tau)) * prod_b (X - j((tau+b)/l)), and each
    X-coefficient -- a level-1 modular function -- is read off as a polynomial in j(tau).

    l must be prime (the l conjugates j((tau+b)/l) are the cyclic l-isogenies)."""
    from nt import primeQ
    if not primeQ(l):
        raise ValueError(f'compute_modpoly needs l prime, got l={l}')
    QMAX = 2*l + 6
    UMAX = l*QMAX + l
    if UMAX > _JQ_NMAX:
        raise ValueError(f'need j-coefficients up to q^{UMAX} for l={l}, have up to q^{_JQ_NMAX} '
                         f'(regenerate data/jq_coeffs.json with more terms)')
    small = {n: _jc(n) for n in range(-1, UMAX+l+2)}      # j(tau/l) in u = q^{1/l}
    # (need terms past UMAX so that small^k stays exact up to u^UMAX for k up to l)
    big   = {l*n: _jc(n) for n in range(-1, QMAX//l + 3)}  # j(l tau), q-series
    jser  = {n: _jc(n) for n in range(-1, QMAX+2)}         # j(tau), q-series
    P = {}
    S = {0: 1}
    for k in range(1, l+1):
        S = _smulser(S, small, UMAX)                      # accumulate j(tau/l)^k incrementally
        Pk = {}
        for e, v in S.items():
            if e % l == 0 and e//l <= QMAX:
                Pk[e//l] = Pk.get(e//l, 0) + l*v
        P[k] = Pk
    sig = {0: {0: 1}}                                      # elementary symmetric of the small conj.
    for k in range(1, l+1):
        acc = {}
        for i in range(1, k+1):
            acc = _sadd(acc, _smul(_smulser(sig[k-i], P[i], QMAX), (-1)**(i-1)))
        sig[k] = {e: v//k for e, v in acc.items() if v % k == 0 and v//k}
    jpows = {0: {0: 1}}
    for d in range(1, l+2):
        jpows[d] = _smulser(jpows[d-1], jser, QMAX)
    def _as_poly_in_j(f):
        f = dict(f); co = {}
        while f and min(f) < 0:
            v = min(f); d = -v; a = f[v]; co[d] = a
            f = _sadd(f, _smul(jpows[d], -a))
        co[0] = f.get(0, 0)
        return co
    M = [[0]*(l+2) for _ in range(l+2)]
    for m in range(0, l+2):                               # coeff of X^{l+1-m} is (-1)^m e_m(j)
        Cm = _smul(_sadd(sig.get(m, {}), _smulser(big, sig.get(m-1, {}), QMAX)), (-1)**m)
        for d, a in _as_poly_in_j(Cm).items():
            M[l+1-m][d] = a
    return M

_MODPOLY_CACHE_PATH = _DATA_DIR / 'classical_modpolys.json'
try:
    with open(_MODPOLY_CACHE_PATH) as f:
        _modpoly_cache = {int(l): M for l, M in json.load(f).items()}
except FileNotFoundError:
    _modpoly_cache = {}

def classical_modpoly(l:int, use_cache=True):
    """Phi_l(X, Y) as a dense (l+2)x(l+2) integer matrix, cached in classical_modpolys.json."""
    if use_cache and l in _modpoly_cache:
        return _modpoly_cache[l]
    M = compute_modpoly(l)
    _modpoly_cache[l] = M
    return M

def save_modpoly_cache(path=_MODPOLY_CACHE_PATH):
    with open(path, 'w') as f:
        json.dump({str(l): M for l, M in _modpoly_cache.items()}, f)

def modpoly_from_terms(l:int, terms):
    """Dense (l+2)x(l+2) Phi_l matrix from a list of (i, j, coeff) terms.  Phi_l is
    symmetric, so each unordered pair may be listed once; both (i,j) and (j,i) are set.
    Use this to ingest a *sourced* classical modular polynomial (e.g. for large l, where
    computing from the q-expansion is slow) into the same cache the computed ones use."""
    M = [[0]*(l+2) for _ in range(l+2)]
    for i, j, co in terms:
        M[i][j] = co
        M[j][i] = co
    return M

def register_modpoly(l:int, M, save=False):
    """Add a Phi_l matrix (computed or sourced) to the in-memory cache; persist if save."""
    _modpoly_cache[l] = M
    if save:
        save_modpoly_cache()

def modpoly_at_j(j:int, l:int, p:int):
    """Coefficients (low-to-high) of Phi_l(j, Y) mod p, a univariate polynomial in Y."""
    M = classical_modpoly(l)
    jp = j % p
    jpows = [1]*(l+2)
    for i in range(1, l+2):
        jpows[i] = (jpows[i-1]*jp) % p
    return [sum(M[i][d]*jpows[i] for i in range(l+2)) % p for d in range(l+2)]

def modpoly_nbrs(j:int, l:int, p:int):
    """The j-invariants l-isogenous to j over F_p: the roots of Phi_l(j, Y) mod p."""
    coeffs = modpoly_at_j(j, l, p)                        # low-to-high
    return [y for y in range(p) if poly_eval_mod(coeffs, y, p, rev=True) == 0]

def modpoly_roots_among(j:int, l:int, p:int, candidates) -> dict:
    """{y: multiplicity} of the roots of Phi_l(j, Y) mod p among candidates.

    Root multiplicity counts the cyclic l-subgroups C with j(E/C) = y, so at
    j = 0 or 1728 it carries the extra-automorphism factor (3 resp. 2) on top
    of the isogeny count."""
    coeffs = modpoly_at_j(j, l, p)                        # low-to-high, deg l+1
    out = {}
    for y in dict.fromkeys(candidates):
        f = coeffs
        m = 0
        while poly_eval_mod(f, y, p, rev=True) == 0:
            m += 1
            q = [0] * (len(f) - 1)                        # deflate by (Y - y)
            acc = 0
            for i in range(len(f) - 1, 0, -1):
                acc = (acc * y + f[i]) % p
                q[i - 1] = acc
            f = q
        if m:
            out[y] = m
    return out

def modpoly_primes() -> list[int]:
    """Primes with a classical Phi_l in the cache (grows via register_modpoly)."""
    return sorted(_modpoly_cache)

def modular_prime_pool() -> list[int]:
    """Every prime with a modular polynomial in SOME format: the 15 genus-0
    Atkin primes plus the classical Phi_l cache.  This is the prime pool the
    j-side neighbour/vertical scans can use, hence the pool for the rigid
    l-set search."""
    return sorted(set(atkin_polys_dict) | set(_modpoly_cache))