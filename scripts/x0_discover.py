"""Discovery-by-evaluation of X_0(ell) plane models and j-maps.

Pipeline (appendix of arXiv:2303.09096, steps 4-6, done numerically):

  1. Sample points (u0, v0) on C_ell(F_p) -- the Tate-normal-form locus where
     P0 = (0,0) has order ell -- by picking random u0 and finding roots of
     r_ell(u0, v) mod p.
  2. At each sample, the diamond operators <m> permute the +-pairs of the
     kernel <P0>, so the power sums  A = sum x([i]P0),  B = sum x([i]P0)^2
     (i = 1..(ell-1)/2) are functions on X_0(ell).  j(E) and the Velu
     codomain j' = j(E/<P0>) = j o Fricke are computed alongside.
  3. Linear algebra over F_p finds the plane relation R(A, B) = 0 and the
     maps j = P(A,B)/Q(A,B), j' = P'(A,B)/Q'(A,B).
  4. Running several primes + CRT + rational reconstruction lifts the
     coefficients to Q; a fresh prime verifies the exact result.

Everything is pure Python; sympy is used only to produce r_ell (x1_derive).
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from x1_derive import eds_sequence, r_ell_fast, u as _u, v as _v  # noqa: E402
import sympy as sp  # noqa: E402


# ---------------------------------------------------------------------------
# r_ell as integer coefficient dicts {(i, j): c}
# ---------------------------------------------------------------------------
_rell_cache = {}


def rell_dict(ell):
    if ell not in _rell_cache:
        W = eds_sequence(ell)
        poly = sp.Poly(r_ell_fast(ell, W), _u, _v)
        _rell_cache[ell] = {m: int(c) for m, c in zip(poly.monoms(), poly.coeffs())}
    return _rell_cache[ell]


# ---------------------------------------------------------------------------
# univariate polynomial arithmetic mod p (dense lists, low -> high)
# ---------------------------------------------------------------------------
def _ptrim(f):
    while f and f[-1] == 0:
        f.pop()
    return f


def _pmulmod(f, g, h, p):
    """f*g mod (h, p); h monic."""
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for k, b in enumerate(g):
                r[i + k] = (r[i + k] + a * b) % p
    return _pdivmod(r, h, p)


def _pdivmod(f, h, p):
    """f mod h (h monic), destructive on a copy."""
    f = f[:]
    dh = len(h) - 1
    while len(f) - 1 >= dh:
        c = f[-1]
        if c:
            sh = len(f) - 1 - dh
            for i, b in enumerate(h):
                f[sh + i] = (f[sh + i] - c * b) % p
        f.pop()
    return _ptrim(f)


def _pgcd(f, g, p):
    f, g = _ptrim(f[:]), _ptrim(g[:])
    while g:
        inv = pow(g[-1], p - 2, p)
        gm = [c * inv % p for c in g]
        f, g = g, _pdivmod(f, gm, p)
    if f:
        inv = pow(f[-1], p - 2, p)
        f = [c * inv % p for c in f]
    return f


def _ppowmod(base, e, h, p):
    r, b = [1], _pdivmod(base, h, p)
    while e:
        if e & 1:
            r = _pmulmod(r, b, h, p)
        b = _pmulmod(b, b, h, p)
        e >>= 1
    return r


def poly_roots_mod(coeffs, p, rng):
    """Roots in F_p of the poly with given low->high coeffs (squarefree-ish)."""
    f = _ptrim([c % p for c in coeffs])
    if not f:
        return []
    inv = pow(f[-1], p - 2, p)
    f = [c * inv % p for c in f]
    # split off the linear factors: gcd(f, x^p - x)
    xp = _ppowmod([0, 1], p, f, p)
    g = xp[:]
    while len(g) < 2:
        g.append(0)
    g[1] = (g[1] - 1) % p
    g = _pgcd(f, _ptrim(g), p)
    roots = []

    def split(h):
        d = len(h) - 1
        if d == 0:
            return
        if d == 1:
            roots.append(-h[0] % p)
            return
        while True:
            a = rng.randrange(p)
            # gcd(h, (x+a)^((p-1)/2) - 1)
            t = _ppowmod([a, 1], (p - 1) // 2, h, p)
            t = t[:]
            if not t:
                t = [0]
            t[0] = (t[0] - 1) % p
            d1 = _pgcd(h, _ptrim(t), p)
            if 0 < len(d1) - 1 < d:
                inv1 = pow(d1[-1], p - 2, p)
                d1 = [c * inv1 % p for c in d1]
                # h / d1
                q = _pquo(h, d1, p)
                split(d1)
                split(q)
                return

    if len(g) > 1:
        split(g)
    return roots


def _pquo(f, h, p):
    """Quotient f // h for monic h, exact division assumed."""
    f = f[:]
    dh = len(h) - 1
    q = [0] * (len(f) - dh)
    while len(f) - 1 >= dh:
        c = f[-1]
        q[len(f) - 1 - dh] = c
        if c:
            sh = len(f) - 1 - dh
            for i, b in enumerate(h):
                f[sh + i] = (f[sh + i] - c * b) % p
        f.pop()
        _ptrim(f)
        if not f:
            break
    return _ptrim(q)


# ---------------------------------------------------------------------------
# long-Weierstrass arithmetic mod p
# ---------------------------------------------------------------------------
class ECp:
    def __init__(self, a1, a2, a3, a4, a6, p):
        self.a = tuple(x % p for x in (a1, a2, a3, a4, a6))
        self.p = p

    def b_invariants(self):
        a1, a2, a3, a4, a6 = self.a
        p = self.p
        b2 = (a1 * a1 + 4 * a2) % p
        b4 = (2 * a4 + a1 * a3) % p
        b6 = (a3 * a3 + 4 * a6) % p
        b8 = (a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4) % p
        return b2, b4, b6, b8

    def disc_j(self):
        p = self.p
        b2, b4, b6, b8 = self.b_invariants()
        disc = (-b2 * b2 % p * b8 - 8 * pow(b4, 3, p) - 27 * b6 * b6 + 9 * b2 * b4 % p * b6) % p
        c4 = (b2 * b2 - 24 * b4) % p
        if disc == 0:
            return 0, None
        return disc, pow(c4, 3, p) * pow(disc, p - 2, p) % p

    def neg(self, P):
        if P is None:
            return None
        a1, a2, a3, a4, a6 = self.a
        x, y = P
        return (x, (-y - a1 * x - a3) % self.p)

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        a1, a2, a3, a4, a6 = self.a
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2 and (y1 + y2 + a1 * x2 + a3) % p == 0:
            return None
        if P == Q:
            lam = (3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1) * \
                pow(2 * y1 + a1 * x1 + a3, p - 2, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
        nu = (y1 - lam * x1) % p
        x3 = (lam * lam + a1 * lam - a2 - x1 - x2) % p
        y3 = (-(lam * x3 + nu) - a1 * x3 - a3) % p
        return (x3, y3)

    def velu_j(self, reps):
        """j-invariant of E/<kernel>, reps = (<P> - O)/{+-1} representatives."""
        a1, a2, a3, a4, a6 = self.a
        p = self.p
        b2 = (a1 * a1 + 4 * a2) % p
        t = w = 0
        for (xq, yq) in reps:
            gx = (3 * xq * xq + 2 * a2 * xq + a4 - a1 * yq) % p
            gy = (-2 * yq - a1 * xq - a3) % p
            tq = (2 * gx - a1 * gy) % p
            uq = gy * gy % p
            t = (t + tq) % p
            w = (w + uq + xq * tq) % p
        A4 = (a4 - 5 * t) % p
        A6 = (a6 - b2 * t - 7 * w) % p
        return ECp(a1, a2, a3, A4, A6, p).disc_j()[1]


def tate_curve(u0, v0, p):
    return ECp(1 - u0, -v0, -v0, 0, 0, p)


def tate_renormalize(E, Q):
    """Tate normal form (u', v') of (E, Q): move Q to (0,0), kill a4, rescale.

    Returns None on degenerate input (Q of low order / vanishing scalings).
    This is the diamond operator when Q = [m]P0."""
    p = E.p
    xQ, yQ = Q
    a1, a2, a3, a4, a6 = E.a
    # shift x -> x + xQ, y -> y + yQ  (r = xQ, s = 0, t = yQ)
    A1 = a1
    A2 = (a2 + 3 * xQ) % p
    A3 = (a3 + xQ * a1 + 2 * yQ) % p
    A4 = (a4 + 2 * xQ * a2 - yQ * a1 + 3 * xQ * xQ) % p
    if A3 == 0:
        return None
    # y -> y + s*x with s = A4/A3 kills a4
    s = A4 * pow(A3, p - 2, p) % p
    A1_, A2_, A3_ = (A1 + 2 * s) % p, (A2 - s * A1 - s * s) % p, A3
    if A2_ == 0:
        return None
    # scale ai -> ai * lam^i with lam = A2_/A3_ to force a2 == a3
    lam = A2_ * pow(A3_, p - 2, p) % p
    u1 = (1 - A1_ * lam) % p
    v1 = (-A2_ * lam * lam) % p
    return u1, v1


# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------
def sample_points(ell, p, count, rng, rd=None):
    """Points (u0, v0) on C_ell(F_p) with smooth fiber, distinct diamond orbits."""
    rd = rd or rell_dict(ell)
    dv = max(m[1] for m in rd)
    out, seen = [], set()
    while len(out) < count:
        u0 = rng.randrange(2, p)
        upow = [1]
        du = max(m[0] for m in rd)
        for _ in range(du):
            upow.append(upow[-1] * u0 % p)
        coeffs = [0] * (dv + 1)
        for (i, jj), c in rd.items():
            coeffs[jj] = (coeffs[jj] + c * upow[i]) % p
        for v0 in poly_roots_mod(coeffs, p, rng):
            E = tate_curve(u0, v0, p)
            disc, jval = E.disc_j()
            if disc == 0:
                continue
            key = (u0, v0)
            if key in seen:
                continue
            seen.add(key)
            out.append((u0, v0))
            if len(out) >= count:
                break
    return out


def invariants(ell, u0, v0, p):
    """(A, B, j, j') at a sample; None if the fiber degenerates."""
    E = tate_curve(u0, v0, p)
    disc, jval = E.disc_j()
    if disc == 0:
        return None
    P0 = (0, 0)
    reps, Q = [], P0
    for _ in range((ell - 1) // 2):
        if Q is None:
            return None          # P0 had smaller order: not on the good locus
        reps.append(Q)
        Q = E.add(Q, P0)
    # Central-moment ratios of the kernel x-multiset: the diamond operators
    # act on {x_i} by an affine map x -> lam^2 (x - xQ), so m_k transforms
    # with weight lam^(2k) and the ratios below are functions on X_0(ell).
    xs = [x for x, _ in reps]
    n = len(xs)
    xbar = sum(xs) * pow(n, p - 2, p) % p
    m2 = sum((x - xbar) ** 2 for x in xs) % p
    m3 = sum((x - xbar) ** 3 for x in xs) % p
    m4 = sum((x - xbar) ** 4 for x in xs) % p
    if m2 == 0 or m3 == 0:
        return None
    A = m4 * pow(m2 * m2 % p, p - 2, p) % p
    B = pow(m2, 3, p) * pow(m3 * m3 % p, p - 2, p) % p
    jp = E.velu_j(reps)
    if jp is None:
        return None
    return A, B, jval, jp


# ---------------------------------------------------------------------------
# linear algebra mod p
# ---------------------------------------------------------------------------
def nullspace_mod(rows, ncols, p):
    """Basis of the right-kernel of the matrix (list of rows) over F_p."""
    m = [r[:] for r in rows]
    nr = len(m)
    pivots = {}
    r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, nr) if m[i][c] % p), None)
        if pr is None:
            continue
        m[r], m[pr] = m[pr], m[r]
        inv = pow(m[r][c], p - 2, p)
        m[r] = [x * inv % p for x in m[r]]
        for i in range(nr):
            if i != r and m[i][c]:
                f = m[i][c]
                m[i] = [(a - f * b) % p for a, b in zip(m[i], m[r])]
        pivots[c] = r
        r += 1
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        vec = [0] * ncols
        vec[fc] = 1
        for c, pr in pivots.items():
            vec[c] = (-m[pr][fc]) % p
        basis.append(vec)
    return basis


def monomials_upto(D, maxb=None):
    return [(a, b) for d in range(D + 1) for a in range(d + 1)
            for b in [d - a] if maxb is None or b <= maxb]


def find_relation(samples, p, maxD=16):
    """Minimal-degree R with R(A,B) = 0 at all samples; returns (D, monos, vec)."""
    for D in range(2, maxD + 1):
        monos = monomials_upto(D)
        if len(samples) < len(monos) + 8:
            raise ValueError('need more samples for degree %d' % D)
        rows = []
        for (A, B, *_rest) in samples:
            Apow = [pow(A, a, p) for a in range(D + 1)]
            Bpow = [pow(B, b, p) for b in range(D + 1)]
            rows.append([Apow[a] * Bpow[b] % p for (a, b) in monos])
        ker = nullspace_mod(rows, len(monos), p)
        if ker:
            return D, monos, ker
    return None


def find_map(samples, vals, p, monos):
    """vec (num_coeffs | den_coeffs) with  N(A,B) - val*Dn(A,B) = 0 at samples.

    monos: list of (a, b) used for both numerator and denominator."""
    rows = []
    for (A, B, *_r), val in zip(samples, vals):
        Apow = {}
        Bpow = {}
        row = []
        for (a, b) in monos:
            Apow.setdefault(a, pow(A, a, p))
            Bpow.setdefault(b, pow(B, b, p))
            row.append(Apow[a] * Bpow[b] % p)
        row += [(-val * x) % p for x in row]
        rows.append(row)
    return nullspace_mod(rows, 2 * len(monos), p)


# ---------------------------------------------------------------------------
# rational reconstruction
# ---------------------------------------------------------------------------
def rat_recon(a, m):
    """r/s = a mod m with |r|, s <= sqrt(m/2); None if it fails."""
    a %= m
    r0, r1 = m, a
    s0, s1 = 0, 1
    bound = int((m // 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if abs(s1) > bound or r1 == 0 and s1 == 0:
        return None
    from math import gcd
    if gcd(r1, s1) != 1:
        return None
    return (r1, s1) if s1 > 0 else (-r1, -s1)


def crt_pair(a1, m1, a2, m2):
    from math import gcd
    assert gcd(m1, m2) == 1
    t = (a2 - a1) * pow(m1, -1, m2) % m2
    return (a1 + m1 * t) % (m1 * m2), m1 * m2


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
PRIMES = [10**9 + 7, 10**9 + 9, 10**9 + 21, 10**9 + 33]


def _discover_one_prime(ell, p, nsamples, seed):
    rng = random.Random(seed)
    pts = sample_points(ell, p, nsamples, rng)
    invs = []
    seen = set()
    for (u0, v0) in pts:
        r = invariants(ell, u0, v0, p)
        if r is None or (r[0], r[1]) in seen:
            continue
        seen.add((r[0], r[1]))
        invs.append(r)
    D, monos, ker = find_relation(invs, p)
    if len(ker) != 1:
        raise RuntimeError('relation kernel dim %d at degree %d' % (len(ker), D))
    rel = ker[0]
    dB = max(b for (a, b), c in zip(monos, rel) if c)
    maps = {}
    for idx, name in ((2, 'j'), (3, 'jprime')):
        vals = [r[idx] for r in invs]
        for Dj in range(2, 3 * D + 4):
            mj = monomials_upto(Dj, maxb=dB - 1)
            if 2 * len(mj) + 8 > len(invs):
                raise RuntimeError('need more samples for %s at degree %d' % (name, Dj))
            kerj = find_map(invs, vals, p, mj)
            if kerj:
                maps[name] = (mj, kerj[0])
                break
        else:
            raise RuntimeError('no %s map found' % name)
    return {'D': D, 'monos': monos, 'rel': rel, 'maps': maps}


def _reconstruct(vec_by_prime, primes):
    """CRT each entry across primes, then rational reconstruction."""
    out = []
    for entries in zip(*vec_by_prime):
        a, m = entries[0] % primes[0], primes[0]
        for e, p in zip(entries[1:], primes[1:]):
            a, m = crt_pair(a, m, e % p, p)
        r = rat_recon(a, m)
        if r is None:
            return None
        out.append(sp.Rational(r[0], r[1]))
    return out


def derive(ell, nsamples=None, nprimes=2, verbose=True):
    """Full discovery for one ell: plane model R(A,B)=0, j and j' maps, exact/QQ.

    Returns dict with sympy expressions in symbols A, B, verified mod a fresh
    prime.  Raises if reconstruction or verification fails (add primes)."""
    nsamples = nsamples or 800
    runs, used = [], []
    for p in PRIMES:
        if len(runs) == nprimes:
            break
        r = _discover_one_prime(ell, p, nsamples, seed=1000 + ell)
        runs.append(r)
        used.append(p)
    r0 = runs[0]
    if any(r['D'] != r0['D'] or len(r['rel']) != len(r0['rel']) for r in runs):
        raise RuntimeError('inconsistent relation shape across primes')
    rel_q = _reconstruct([r['rel'] for r in runs], used)
    maps_q = {}
    for name in ('j', 'jprime'):
        if any(len(r['maps'][name][1]) != len(r0['maps'][name][1]) for r in runs):
            raise RuntimeError('inconsistent %s map shape' % name)
        maps_q[name] = _reconstruct([r['maps'][name][1] for r in runs], used)
    if rel_q is None or any(m is None for m in maps_q.values()):
        raise RuntimeError('rational reconstruction failed; increase nprimes')

    A, B = sp.symbols('A B')

    def poly_of(monos, coeffs):
        return sp.Add(*[c * A**a * B**b for (a, b), c in zip(monos, coeffs) if c])

    R = poly_of(r0['monos'], rel_q)
    # clear denominators, make primitive
    R = sp.Poly(R * sp.lcm([sp.fraction(c)[1] for c in rel_q if c]), A, B).primitive()[1].as_expr()
    result = {'ell': ell, 'R': R}
    for name in ('j', 'jprime'):
        mj = r0['maps'][name][0]
        cq = maps_q[name]
        n = len(mj)
        num, den = poly_of(mj, cq[:n]), poly_of(mj, [-c for c in cq[n:]])
        g = sp.gcd(sp.Poly(num, A, B), sp.Poly(den, A, B))
        result[name] = sp.cancel((num / g.as_expr()) / (den / g.as_expr()))
    ok = verify(result, nsamples=40)
    result['verified'] = ok
    if verbose:
        print('ell=%d: relation degree %d, verified %s' % (ell, r0['D'], ok))
    return result


def verify(model, nsamples=40, p=None):
    """Fresh-prime check: R(A,B)=0, j and j' match Velu ground truth."""
    ell = model['ell']
    p = p or 999999937
    rng = random.Random(31337)
    A, B = sp.symbols('A B')
    Rp = sp.Poly(model['R'], A, B)
    pts = sample_points(ell, p, nsamples, rng)
    checked = 0
    for (u0, v0) in pts:
        r = invariants(ell, u0, v0, p)
        if r is None:
            continue
        Av, Bv, jv, jpv = r
        if int(Rp.eval((Av, Bv))) % p != 0:
            return False
        for name, truth in (('j', jv), ('jprime', jpv)):
            num, den = sp.fraction(model[name])
            nv = int(sp.Poly(num, A, B).eval((Av, Bv))) % p
            dv = int(sp.Poly(den, A, B).eval((Av, Bv))) % p
            if dv == 0:
                continue
            if nv * pow(dv, p - 2, p) % p != truth:
                return False
        checked += 1
    return checked >= nsamples // 2


if __name__ == '__main__':
    for ell in [int(x) for x in sys.argv[1:]] or [17, 19]:
        m = derive(ell)
        print('R_%d(A,B) =' % ell, m['R'])
        print('j_%d =' % ell, m['j'])
