"""Elliptic curves over F_p, by F_p methods only.

Model/j-invariant conversions, the signature invariant, fast trace-of-Frobenius
tables, and Atkin-polynomial isogeny codomains.  Everything here works purely
over the prime field -- no lattice or quadratic-form data.

Planned direction: this module is the home for a proper elliptic-curve-over-F_p
class (point arithmetic, group structure) and eventually the Velu isogeny
engine now in velu.py.

History note: the pre-ecqf_bij bijection path (sig_to_qfs_dict,
get_j_to_qfs_dict, the frobmat/mw_gens kernel chain, and the endomorphism-disc
candidate search) lived here until the 2026-08 cleanup; the working successors
are ecqf_bij.ecqf_full_bijection_ord and the ecqf_tools Mordell-Weil chain
(qf_ap_FrMat / frob_to_mw_gens).  See tag `pre-cleanup` for the old code.
"""

from nt import quad_rec, discfac, find_nonsquare

import numpy as np
from functools import lru_cache


def j_to_fg(j:int,char = 0):
    if j == 0:
        return (0,1)
    elif j == 1728 or char>0 and (j-1728)%char == 0:
        return (1,0)
    else:
        f = -3*j*(j-1728)
        g = 2*j*((j-1728)**2)
        if char == 0:
            return (f,g)
        else:
            return (f % char, g% char)

def fg_to_j(fg:tuple[int,int],char =0):
    f,g = fg
    if f == 0 or char>0 and f%char ==0:
        return 0
    elif g == 0 or char>0 and g % char ==0:
        return 1728
    else:
        f3 = 4*(f**3)
        jnum = 1728 *f3
        jden = f3+27*(g**2)
        if char == 0:
            if jden == 0:
                raise ZeroDivisionError('Singular curve')
            elif jnum % jden == 0:
                return jnum//jden
            else:
                return jnum/jden
        else:
            jnum = jnum % char
            jden = jden % char
            if jden == 0:
                raise ZeroDivisionError('Singular curve')
            jdeninv = pow(jden,-1,char)
            return (jnum*jdeninv)%char

def signature(f:int,g:int,p:int)->int:
    """Signature s in {+1,-1} of the model y^2 = x^3 + f x + g over F_p -- the
    F_p-isomorphism invariant that distinguishes a supersingular curve from its
    quadratic twist.  For j != 1728 it is the QR character of the constant term g;
    for j == 1728 (g == 0) the QR character of the x-coefficient f.  (Under
    (f,g) |-> (u^4 f, u^6 g) both u^4 and u^6 are squares, so the character is an
    iso invariant, and a non-square twist u flips it.)"""
    f %= p; g %= p
    return quad_rec(f,p) if g == 0 else quad_rec(g,p)

def js_to_fg(js:tuple[int,int],p:int)->tuple[int,int]:
    """Canonical short Weierstrass model (f,g) for the signature (j, s) over F_p,
    i.e. a model with j-invariant j (mod p) and signature(f,g,p) == s.  Mirrors
    j_to_fg but selects the curve-vs-twist representative named by s, with the
    fixed conventions y^2 = x^3 + s*x at j == 1728 and (0,1)/(0,-3) at j == 0."""
    j,s = js
    if (j-1728) % p == 0:
        return (s % p, 0)
    if j % p == 0:
        return (0,1) if s == 1 else (0,(-3) % p)
    f,g = j_to_fg(j % p, p)
    if quad_rec(g,p) == s:
        return (f % p, g % p)
    t = find_nonsquare(p)
    return ((pow(t,2,p)*f) % p, (pow(t,3,p)*g) % p)

@lru_cache(maxsize=8)
def _chi_table(p:int):
    """Quadratic character of F_p as a length-p int8 lookup table (chi[0] = 0)."""
    chi = np.full(p, -1, dtype=np.int8)
    x = np.arange(p, dtype=np.int64)
    chi[(x * x) % p] = 1
    chi[0] = 0
    return chi


def trace_frob(fg:tuple[int,int],p:int)->int:
    f,g = fg
    if p < 5:                                 # tiny p: keep quad_rec's conventions
        return - sum([quad_rec(x**3+f*x+g,p) for x in range(p)])
    x = np.arange(p, dtype=np.int64)
    vals = ((x * x % p) * x + f * x + g) % p  # x^3+fx+g, staying inside int64
    return - int(_chi_table(p)[vals].sum())


def fp_isog_codomains(j:int,l:int,p:int):
    """j-invariants l-isogenous to j over F_p, via the Atkin modular polynomial.

    Roots of the Atkin polynomial at j are found by a vectorized Horner scan
    over all of F_p; each root y0 maps to the codomain a(y0) - j.  (The
    modularpolynomials import is deferred so importing ecfp stays cheap --
    that module loads its polynomial stores eagerly.)"""
    from modularpolynomials import atk_at_j_coeffs, atk_a_coeffs, poly_eval_mod
    coeffs = atk_at_j_coeffs(j, l, p)
    x = np.arange(p, dtype=np.int64)
    vals = np.zeros(p, dtype=np.int64)
    for c in reversed(coeffs):
        vals = (vals * x + c) % p
    ca = atk_a_coeffs(l)
    return [(poly_eval_mod(ca[::-1], int(y0), p) - j) % p
            for y0 in np.flatnonzero(vals == 0)]


TRACE_TABLE_IMPL = None   # optional faster backend; see trace_gpu.enable()


@lru_cache(maxsize=8)
def _trace_table(p:int):
    """tr[j] = trace of Frobenius of the canonical model j_to_fg(j) over F_p, for
    every j at once (blocked numpy scan).  Entries at j = 0 and j = 1728 are for
    the singular formula output and must not be used -- the callers special-case
    those two j's.  Cached per p: every class at p, and every re-scan within one
    bijection computation, reads the same table.  A plugged-in backend
    (TRACE_TABLE_IMPL, e.g. the GPU kernel in trace_gpu) may take over; returning
    None from it falls back to the numpy scan."""
    if TRACE_TABLE_IMPL is not None:
        tr = TRACE_TABLE_IMPL(p)
        if tr is not None:
            tr = np.asarray(tr, dtype=np.int64)
            tr.setflags(write=False)
            return tr
    x = np.arange(p, dtype=np.int64)
    cubes = (x * x % p) * x % p
    chi = _chi_table(p)
    t = (x - 1728) % p
    f = (-3 * x % p) * t % p                               # j_to_fg, vectorized
    g = (2 * x % p) * (t * t % p) % p
    tr = np.empty(p, dtype=np.int64)
    B = max(1, (1 << 22) // p)
    for s in range(0, p, B):
        e = min(s + B, p)
        vals = (cubes[None, :] + f[s:e, None] * x[None, :] + g[s:e, None]) % p
        tr[s:e] = -chi[vals].sum(axis=1)
    tr.setflags(write=False)
    return tr


def trfr_to_js(a:int,p:int):
    if p < 5:
        js = [j for j in range(1,p) if
                (j-1728)%p!=0 and abs(trace_frob(j_to_fg(j),p))==abs(a)]
    else:
        tr = _trace_table(p)
        js = [j for j in range(1, p) if (j - 1728) % p != 0 and abs(int(tr[j])) == abs(a)]
    d = discfac(a*a-4*p)[0]
    if d==-3 or (d%p == 0 and p % 3 == 2):
        js.append(0)
    if d == -4 or (d % p == 0 and p % 4 == 3):
        js.append(1728%p)
    return js

def trfr_to_models(a:int,p:int):
    d,c = discfac(a**2-4*p)
    j_to_fg_dict = {}
    models = []
    # We need a nonsquare mod p
    ns = p-1
    while ns>1 and quad_rec(ns,p)>0:
        ns-=1
    if d == -3:
        fg0 = (0,1)
        tr0 = trace_frob(fg0,p)
        if tr0 == a:
            j_to_fg_dict[0]=fg0
        elif tr0 == -a:
            j_to_fg_dict[0] = (0,pow(ns,3,p))
        else:
            gn = p-2
            while pow(gn,p//3,p)==1:
                gn-=1
            tr1 = trace_frob((0,gn),p)
            if tr1 == a:
                j_to_fg_dict[0] = (0,gn)
            elif tr1 == -a:
                j_to_fg_dict[0] = (0,(pow(ns,3,p)*gn)%p)
            else:
                tr2 = trace_frob((0,pow(gn,2,p)),p)
                if tr2 == a:
                    j_to_fg_dict[0]=(0,pow(gn,2,p))
                else:
                    j_to_fg_dict[0] = (0,(pow(ns,3,p)*pow(gn,2,p))%p)
    elif d == -4:
        fg0 = (1,0)
        tr0 = trace_frob(fg0,p)
        if tr0 == a:
            j_to_fg_dict[1728%p]=fg0
        elif tr0 == -a:
            j_to_fg_dict[1728%p]=(pow(ns,2,p),p)
        else:
            tr1 = trace_frob((ns,0),p)
            if tr1 == a:
                j_to_fg_dict[1728%p] = (ns,0)
            else:
                j_to_fg_dict[1728%p] = (pow(ns,3,p),0)
    elif a %p == 0:
        if p % 3 == 2:
            models+=[(0,1),(0,ns)]
        if p % 4 == 3:
            models+=[(p-1,0),(1,0)]
    table = _trace_table(p) if p >= 5 else None
    for j0 in range(1,p):
        if (j0-1728)%p!=0:
            fg0 = j_to_fg(j0,p)
            tr0 = int(table[j0]) if table is not None else trace_frob(fg0,p)
            if tr0 == a:
                if a == 0:
                    f0, g0 = fg0
                    j_to_fg_dict = [fg0,((f0*(ns**2))%p,(g0*(ns**3))%p)]
                else:
                    j_to_fg_dict[j0] = fg0
                models.append(fg0)
            elif tr0 == -a:
                f0,g0 = fg0
                j_to_fg_dict[j0]=((f0*(ns**2))%p,(g0*(ns**3))%p)
    return j_to_fg_dict
