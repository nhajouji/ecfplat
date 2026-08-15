"""Discover the quotient map C(ell) -> X_0(ell)-model (Nadir's curves) by evaluation.

Per sample (u,v) on C(ell)(F_p): compute (j, j') from the universal curve; find the
unique (x,y) on the target model with j_model(x,y) = j and j_model(Fricke) = j';
then solve for x(u,v), y(u,v) as rational functions mod r_ell by nullspace
linear algebra, escalating primes until rational reconstruction verifies.
"""
import os, pickle, random, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import x0_discover as X
import sympy as sp

x, y, jsym = sp.symbols('x y j')

MODELS = {
    17: {
        'a': (1, -1, 1, -1, -14),
        'jnum': (24 + 12*x + 6*x**2 - 5*x**3 - 4*y + 8*x*y + x**2*y)**3,
        'jden': -43 + 14*x + x**2 - 8*y,
        'S': (7, 13),          # Fricke: P -> S - P
    },
    19: {
        'a': (0, 1, 1, -9, -15),
        'jnum': (-14 - x)*(-3694 + 2411*x + 1936*x**2 + 96*x**3 + 3840*y + 600*x*y)**3,
        'jden': (x - 5)*(4*x - 11 - y)**6,
        'S': None,             # Fricke: P -> -P
        # O is not a cusp on this model.  Discover in the chart translated by
        # the width-19 cusp c: g(P), h(P) = coords of P - c, whose poles sit
        # only over the cusp family at u in {0, inf} upstairs -> monomial
        # denominators, like ell=17.  Final map: (x,y) = (g,h) + c.
        'chart': ('translate', (5, 9)),
    },
}


def reduced_condition(ell):
    """jnum - j*jden reduced mod the curve to A0(x,j) + A1(x,j)*y; returns
    (A0, A1) as sympy Polys in (x, j) with integer coefficients."""
    cfg = MODELS[ell]
    a1, a2, a3, a4, a6 = cfg['a']
    y2 = x**3 + a2*x**2 + a4*x + a6 - a1*x*y - a3*y   # y^2 == this, linear in y
    expr = sp.expand(cfg['jnum'] - jsym*cfg['jden'])
    P = sp.Poly(expr, y)
    while P.degree() > 1:
        cs = P.all_coeffs()          # high -> low in y
        d = P.degree()
        lead = cs[0]
        rest = sp.Add(*[c*y**(d-1-i) for i, c in enumerate(cs[1:])])
        P = sp.Poly(sp.expand(rest + lead*y**(d-2)*y2), y)
    A1 = sp.Poly(P.coeff_monomial(y) if P.degree() >= 1 else 0, x, jsym)
    A0 = sp.Poly(P.coeff_monomial(1), x, jsym)
    return A0, A1


def _eval_xpoly(P, jval, p):
    """P(x, j) with j = jval, as dense mod-p coeff list (low->high in x)."""
    d = P.degree(x)
    out = [0]*(d+1)
    for (i, k), c in zip(P.monoms(), P.coeffs()):
        out[i] = (out[i] + int(c)*pow(jval, k, p)) % p
    return out


def _pmul(f, g, p):
    r = [0]*(len(f)+len(g)-1)
    for i, a in enumerate(f):
        if a:
            for k, b in enumerate(g):
                r[i+k] = (r[i+k] + a*b) % p
    return r


def _padd(f, g, p):
    n = max(len(f), len(g))
    return [((f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % p
            for i in range(n)]


def match_point(ell, jval, jpval, A0P, A1P, E, p, rng):
    """The point (x0,y0) on the model with j=jval and j(Fricke)=jpval, or None."""
    cfg = MODELS[ell]
    a1, a2, a3, a4, a6 = cfg['a']
    A0 = _eval_xpoly(A0P, jval, p)
    A1 = _eval_xpoly(A1P, jval, p)
    # curve with y = -A0/A1, cleared: A0^2 - a1 x A0 A1 - a3 A0 A1 - (x^3+a2x^2+a4x+a6) A1^2
    F = _pmul(A0, A0, p)
    t = _pmul(A0, A1, p)
    F = _padd(F, [(-a1*c) % p for c in ([0]+t)], p)      # -a1*x*A0*A1
    F = _padd(F, [(-a3*c) % p for c in t], p)
    cube = [(a6) % p, a4 % p, a2 % p, 1]
    F = _padd(F, [(-c) % p for c in _pmul(cube, _pmul(A1, A1, p), p)], p)
    cands = []
    for x0 in X.poly_roots_mod(F, p, rng):
        d = 0
        for i, c in enumerate(A1):
            d = (d + c*pow(x0, i, p)) % p
        if d == 0:
            continue
        n = 0
        for i, c in enumerate(A0):
            n = (n + c*pow(x0, i, p)) % p
        y0 = (-n) * pow(d, p-2, p) % p
        P = (x0, y0)
        if not E.on(P):
            continue
        Q = E.add(cfg['S'], E.neg(P)) if cfg['S'] else E.neg(P)
        jq = model_j(ell, Q, p)
        if jq == jpval:
            cands.append(P)
    return cands[0] if len(cands) == 1 else None


_modelj_cache = {}


def model_j(ell, P, p):
    if ell not in _modelj_cache:
        cfg = MODELS[ell]
        _modelj_cache[ell] = (sp.Poly(sp.expand(cfg['jnum']), x, y),
                              sp.Poly(sp.expand(cfg['jden']), x, y))
    NP, DP = _modelj_cache[ell]
    if P is None:
        return None
    x0, y0 = P

    def ev(Q):
        s = 0
        for (i, k), c in zip(Q.monoms(), Q.coeffs()):
            s = (s + int(c)*pow(x0, i, p)*pow(y0, k, p)) % p
        return s
    d = ev(DP)
    if d == 0:
        return None
    return ev(NP) * pow(d, p-2, p) % p


def matched_samples(ell, p, count, seed, A0P, A1P):
    rng = random.Random(seed)
    cfg = MODELS[ell]
    E = X.ECp(*cfg['a'], p)
    pts = X.sample_points(ell, p, count, rng)
    out = []
    for (u0, v0) in pts:
        r = X.invariants(ell, u0, v0, p)
        if r is None:
            continue
        _, _, jval, jpval = r
        P = match_point(ell, jval, jpval, A0P, A1P, E, p, rng)
        if P is None:
            continue
        x0, y0 = P
        ch = cfg.get('chart')
        if ch and ch[0] == 'translate':
            Q = E.add((x0, y0), E.neg(ch[1]))
            if Q is None:
                continue
            x0, y0 = Q
        out.append((u0, v0, x0, y0))
    return out


def canonical_reduce(vec, mn, md, p):
    """Kernel vector -> primitive canonical rep: strip gcd_u, make d(u) monic.

    Returns ({(a,b): c} numerator, [d0..dk] denominator), canonical across
    primes whatever degree overshoot the search used."""
    n = len(mn)
    numd = {ab: c for ab, c in zip(mn, vec[:n]) if c}
    den = [0]*(1 + max(a for a, _ in md))
    for (a, _), c in zip(md, vec[n:]):
        den[a] = c
    den = X._ptrim(den)
    bmax = max(b for _, b in numd) if numd else 0
    cols = []
    for b in range(bmax + 1):
        amax = max((a for (a, bb) in numd if bb == b), default=-1)
        cols.append([numd.get((a, b), 0) for a in range(amax + 1)])
    g = den[:]
    for c in cols:
        if c:
            g = X._pgcd(g, c, p)
    if len(g) > 1:
        den = X._pquo(den, g, p)
        cols = [X._pquo(c, g, p) if c else c for c in cols]
    lead = pow(den[-1], p - 2, p)
    den = [x*lead % p for x in den]
    out = {}
    for b, c in enumerate(cols):
        for a, x in enumerate(c):
            if x:
                out[(a, b)] = x*lead % p
    return out, den


def discover_map(ell, p, samples, dv, Da0=12, maxDa=40):
    """x(u,v), y(u,v) as canonical N(u,v)/d(u) (v-deg < dv in N)."""
    res = {}
    for name, idx in (('x', 2), ('y', 3)):
        vals = [s[idx] for s in samples]
        found = None
        Da = Da0
        while Da <= maxDa:
            mn = [(a, b) for b in range(dv) for a in range(Da + 1)]
            md = [(a, 0) for a in range(Da + 1)]
            if len(mn) + len(md) + 8 > len(samples):
                raise RuntimeError('need more samples at Da=%d' % Da)
            ker = X.find_map(samples, vals, p, mn, md)
            if ker:
                found = canonical_reduce(ker[0], mn, md, p)
                break
            Da += 4
        if not found:
            raise RuntimeError('no %s-map found for ell=%d (maxDa=%d)' % (name, ell, maxDa))
        res[name] = found
    return res


def run(ell, nsamples=1400):
    t0 = time.time()
    cachedir = Path(os.environ.get('QMAP_CACHE',
                                   os.path.expanduser('~/.cache/ecfplat_qmap')))
    cachedir.mkdir(parents=True, exist_ok=True)
    A0P, A1P = reduced_condition(ell)
    rd = X.rell_dict(ell)
    dv = max(m[1] for m in rd)
    runs, used = [], []
    for p in X.PRIMES:
        f = cachedir / f'qmap_{ell}_{p}.pkl'
        if f.exists():
            r = pickle.load(open(f, 'rb'))
        else:
            samples = matched_samples(ell, p, nsamples, 5000 + ell, A0P, A1P)
            r = discover_map(ell, p, samples, dv)
            r['nsamples'] = len(samples)
            pickle.dump(r, open(f, 'wb'))
        shape = {n: (sorted(r[n][0]) and (max(a for a, _ in r[n][0]),
                     max(b for _, b in r[n][0]), len(r[n][1]))) for n in ('x', 'y')}
        print(f'p={p}: {r["nsamples"]} matched, x shape {shape["x"]}, y shape {shape["y"]} '
              f'({time.time()-t0:.0f}s)', flush=True)
        runs.append(r)
        used.append(p)
        if len(runs) < 2:
            continue
        ok_shape = all(set(rr[n][0]) == set(runs[0][n][0]) and len(rr[n][1]) == len(runs[0][n][1])
                       for rr in runs for n in ('x', 'y'))
        if not ok_shape:
            print('  support mismatch across primes; continuing', flush=True)
            continue
        maps_q = {}
        for name in ('x', 'y'):
            keys = sorted(runs[0][name][0])
            vecs = [[rr[name][0][k] for k in keys] + rr[name][1] for rr in runs]
            maps_q[name] = (keys, X._reconstruct(vecs, used))
        if any(v[1] is None for v in maps_q.values()):
            print('  not yet reconstructible', flush=True)
            continue
        print(f'reconstructed with {len(runs)} primes', flush=True)
        u, v = sp.symbols('u v')
        out = {'ell': ell}
        for name in ('x', 'y'):
            keys, cq = maps_q[name]
            n = len(keys)
            num = sp.Add(*[c*u**a*v**b for (a, b), c in zip(keys, cq[:n]) if c])
            den = sp.Add(*[c*u**a for a, c in enumerate(cq[n:]) if c])
            out[name] = num/den
        ch = MODELS[ell].get('chart')
        if ch and ch[0] == 'translate':
            G, H = out['x'], out['y']
            a1, a2, a3, a4, a6 = MODELS[ell]['a']
            cx, cy = ch[1]
            lam = sp.cancel((cy - H)/(cx - G))
            nu = sp.cancel(H - lam*G)
            x3 = sp.cancel(lam**2 + a1*lam - a2 - G - cx)
            out['x'] = x3
            out['y'] = sp.cancel(-(lam*x3 + nu) - a1*x3 - a3)
        ok = verify_map(ell, out, A0P, A1P)
        print('fresh-prime verification:', ok, flush=True)
        if ok:
            out['verified'] = True
            pickle.dump(out, open(cachedir / f'qmap_{ell}.pkl', 'wb'))
            return out
    raise RuntimeError('ran out of primes')


def verify_map(ell, out, A0P, A1P, p=999999937, n=60):
    rng = random.Random(424242)
    cfg = MODELS[ell]
    E = X.ECp(*cfg['a'], p)
    u, v = sp.symbols('u v')
    polys = {}
    for name in ('x', 'y'):
        nm, dn = sp.fraction(sp.together(out[name]))
        polys[name] = (sp.Poly(nm, u, v), sp.Poly(dn, u, v))
    pts = X.sample_points(ell, p, n, rng)
    good = 0
    for (u0, v0) in pts:
        r = X.invariants(ell, u0, v0, p)
        if r is None:
            continue
        _, _, jval, jpval = r
        vals = {}
        skip = False
        for name in ('x', 'y'):
            nm, dn = polys[name]
            dv_ = int(dn.eval((u0, v0))) % p
            if dv_ == 0:
                skip = True
                break
            vals[name] = int(nm.eval((u0, v0))) * pow(dv_, p-2, p) % p
        if skip:
            continue
        P = (vals['x'], vals['y'])
        if not E.on(P):
            return False
        if model_j(ell, P, p) != jval:
            return False
        Q = E.add(cfg['S'], E.neg(P)) if cfg['S'] else E.neg(P)
        if model_j(ell, Q, p) != jpval:
            return False
        good += 1
    return good >= n // 3


if __name__ == '__main__':
    for ell in [int(a) for a in sys.argv[1:]] or [17]:
        out = run(ell)
        u, v = sp.symbols('u v')
        for name in ('x', 'y'):
            nm, dn = sp.fraction(sp.together(out[name]))
            print(f'{name}(u,v): num {len(sp.Poly(nm,u,v).terms())} terms, '
                  f'den {len(sp.Poly(dn,u,v).terms())} terms, '
                  f'max coeff digits {max(len(str(abs(c))) for c in sp.Poly(nm,u,v).coeffs())}')
