# p <-> l duality experiments: symmetry of quadratic-factor counts on the
# supersingular locus of X_0(l)^+.
#
# Stage 1 replicates the nakaya.nb Mathematica computation: for primes p, l in
# the 15 genus-0 (Monster) primes, let S_p = supersingular j-invariants in F_p.
# With (A_l, B_l) the Atkin modular relation j^2 - A_l(t) j + B_l(t) = 0 on
# X_0(l)^+, form
#     G_{p,l}(t) = rad( prod_{j in S_p} (j^2 - A_l(t) j + B_l(t)) mod p )
# and count irreducible factors of G by degree.  Q(p,l) = # quadratic factors.
# Observed (nakaya.nb): Q(p,l) = Q(l,p) for all pairs of Monster primes.
#
# Polynomials over F_p are numpy int64 coefficient arrays, low-to-high degree.

import json
import time
from pathlib import Path

import numpy as np

_DATA_DIR = Path(__file__).parent / 'data'

MONSTER_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]

# Q(p,l) from nakaya.nb (Out[92]): NB_QMAT[r][c] = Q(p = MONSTER_PRIMES[c],
# l = MONSTER_PRIMES[r]).  Symmetric; diagonal is degenerate (set to 0 there).
NB_QMAT = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 1, 1, 2, 2],
    [0, 0, 0, 0, 0, 1, 2, 0, 2, 1, 1, 2, 5, 3, 4],
    [0, 0, 0, 0, 1, 1, 0, 2, 2, 3, 3, 2, 3, 5, 9],
    [0, 0, 0, 1, 0, 1, 3, 2, 5, 4, 4, 6, 10, 12, 15],
    [0, 0, 1, 1, 1, 0, 2, 4, 2, 5, 8, 11, 9, 12, 16],
    [0, 0, 2, 0, 3, 2, 0, 5, 5, 10, 7, 15, 12, 18, 22],
    [0, 1, 0, 2, 2, 4, 5, 0, 7, 5, 10, 11, 16, 17, 20],
    [0, 0, 2, 2, 5, 2, 5, 7, 0, 12, 12, 16, 21, 27, 33],
    [1, 0, 1, 3, 4, 5, 10, 5, 12, 0, 12, 24, 26, 28, 38],
    [0, 2, 1, 3, 4, 8, 7, 10, 12, 12, 0, 17, 28, 34, 39],
    [1, 1, 2, 2, 6, 11, 15, 11, 16, 24, 17, 0, 36, 46, 51],
    [0, 1, 5, 3, 10, 9, 12, 16, 21, 26, 28, 36, 0, 56, 67],
    [1, 2, 3, 5, 12, 12, 18, 17, 27, 28, 34, 46, 56, 0, 83],
    [1, 2, 4, 9, 15, 16, 22, 20, 33, 38, 39, 51, 67, 83, 0],
]

###########################
# F_p polynomial toolkit  #
###########################
# A polynomial is a 1-d np.int64 array of coefficients low-to-high, reduced
# mod p, with no trailing zeros (the zero polynomial is the empty array).

def pnorm(a, p):
    a = np.asarray(a, dtype=np.int64) % p
    nz = np.nonzero(a)[0]
    return a[: nz[-1] + 1] if nz.size else a[:0]

def pdeg(a):
    return len(a) - 1        # zero polynomial gets degree -1

def pmul(a, b, p):
    if len(a) == 0 or len(b) == 0:
        return a[:0]
    # coefficients < p, lengths bounded so int64 convolution is exact
    return np.convolve(a, b) % p

def pmonic(a, p):
    if len(a) == 0:
        return a
    inv = pow(int(a[-1]), p - 2, p)
    return (a * inv) % p

def pdivmod(a, b, p):
    """divmod of a by b (b nonzero), over F_p."""
    if len(b) == 0:
        raise ZeroDivisionError
    a = a.copy()
    db, da = pdeg(b), pdeg(a)
    if da < db:
        return a[:0], a
    inv = pow(int(b[-1]), p - 2, p)
    q = np.zeros(da - db + 1, dtype=np.int64)
    for i in range(da - db, -1, -1):
        c = (a[i + db] * inv) % p
        if c:
            q[i] = c
            a[i : i + db + 1] = (a[i : i + db + 1] - c * b) % p
    return pnorm(q, p), pnorm(a, p)

def pmod(a, b, p):
    return pdivmod(a, b, p)[1]

def pgcd(a, b, p):
    while len(b):
        a, b = b, pmod(a, b, p)
    return pmonic(a, p)

def pmulmod(a, b, f, p):
    return pmod(pmul(a, b, p), f, p)

def ppowmod(a, e, f, p):
    """a^e mod (f, p) by binary exponentiation."""
    result = pnorm(np.array([1]), p)
    a = pmod(a, f, p)
    while e:
        if e & 1:
            result = pmulmod(result, a, f, p)
        a = pmulmod(a, a, f, p)
        e >>= 1
    return result

def pderiv(a, p):
    if len(a) <= 1:
        return a[:0]
    return pnorm(a[1:] * np.arange(1, len(a), dtype=np.int64), p)

def pradical(w, p):
    """Radical (squarefree part) of w, correct in characteristic p.

    Factors f^e with p | e hide entirely inside gcd(w, w'); they are recovered
    from the p-th-power part via the F_p identity h(x)^p = h(x^p)."""
    w = pmonic(w, p)
    if pdeg(w) <= 0:
        return w
    dw = pderiv(w, p)
    if len(dw) == 0:
        return pradical(pnorm(w[::p], p), p)     # w is a p-th power
    g = pgcd(w, dw, p)
    w1 = pdivmod(w, g, p)[0]        # rad of the part with exponent prime to p
    # strip w1's factors out of g; what survives is a p-th power
    while True:
        c = pgcd(g, w1, p)
        if pdeg(c) == 0:
            break
        g = pdivmod(g, c, p)[0]
    if pdeg(g) > 0:
        w1 = pmul(w1, pradical(pnorm(g[::p], p), p), p)
    return pmonic(w1, p)

def ddf_counts(g, p, maxdeg=None):
    """Distinct-degree census of squarefree monic g over F_p.

    Returns {d: number of irreducible factors of degree d}."""
    g = pmonic(g, p)
    counts = {}
    x = pnorm(np.array([0, 1]), p)
    h = x.copy()                      # h = x^(p^d) mod g
    d = 0
    while pdeg(g) > 0:
        d += 1
        if maxdeg is not None and d > maxdeg:
            break
        if 2 * d > pdeg(g):
            # remaining g is irreducible
            counts[pdeg(g)] = counts.get(pdeg(g), 0) + 1
            break
        h = ppowmod(h, p, g, p)
        gd = pgcd(g, psub(h, x, p), p)
        if pdeg(gd) > 0:
            counts[d] = pdeg(gd) // d
            g = pdivmod(g, gd, p)[0]
            h = pmod(h, g, p)
    return counts

def psub(a, b, p):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=np.int64)
    out[: len(a)] += a
    out[: len(b)] -= b
    return pnorm(out, p)

##############################
# Stage 1: notebook replica  #
##############################

_ss_sigs = None

def _load_ss_sigs():
    global _ss_sigs
    if _ss_sigs is None:
        with open(_DATA_DIR / 'ecqf_ss_pcbij_velu_4_1024.json') as f:
            raw = json.load(f)
        _ss_sigs = {int(ps): sorted({int(k[1:-1].split(',')[0]) for k in raw[ps]})
                    for ps in raw}
    return _ss_sigs

def rational_ss_js(p):
    """Supersingular j-invariants lying in F_p."""
    if p == 2 or p == 3:
        return [0]
    sigs = _load_ss_sigs()
    if p not in sigs:
        raise ValueError(f'no precomputed supersingular data for p = {p}')
    return sigs[p]

_atkin = None

def _load_atkin():
    global _atkin
    if _atkin is None:
        with open(_DATA_DIR / 'atkinpolys.json') as f:
            raw = json.load(f)
        _atkin = {int(l): v for l, v in raw.items()}
    return _atkin

def atkin_AB(l, p):
    """Coefficient arrays of A_l, B_l = B1 * B3^3 mod p (low-to-high)."""
    atk = _load_atkin()
    if l not in atk:
        raise ValueError(f'no Atkin polynomial for l = {l}')
    a_c, b1_c, b3_c = atk[l]
    A = pnorm(np.array([c % p for c in a_c]), p)
    B1 = pnorm(np.array([c % p for c in b1_c]), p)
    B3 = pnorm(np.array([c % p for c in b3_c]), p)
    B = pmul(B1, pmul(B3, pmul(B3, B3, p), p), p)
    return A, B

def sspl_radical(p, l):
    """rad( prod_{j in S_p^rat} (j^2 - A_l(t) j + B_l(t)) mod p )."""
    A, B = atkin_AB(l, p)
    G = pnorm(np.array([1]), p)
    for j in rational_ss_js(p):
        j = j % p
        term = psub(B, pnorm(A * j, p), p)
        term = padd(term, pnorm(np.array([(j * j) % p]), p), p)
        G = pmul(G, term, p)
    return pradical(G, p)

def padd(a, b, p):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=np.int64)
    out[: len(a)] += a
    out[: len(b)] += b
    return pnorm(out, p)

def q_profile(p, l, maxdeg=None):
    """Degree census of the irreducible factors of sspl_radical(p, l)."""
    return ddf_counts(sspl_radical(p, l), p, maxdeg=maxdeg)

def Q(p, l):
    """Number of irreducible quadratic factors (the nakaya.nb count)."""
    return q_profile(p, l, maxdeg=2).get(2, 0)

def stage1(verbose=True):
    """Recompute the 15x15 matrix and compare with the notebook's."""
    n = len(MONSTER_PRIMES)
    mat = [[0] * n for _ in range(n)]
    lin = [[0] * n for _ in range(n)]
    t0 = time.time()
    for r, l in enumerate(MONSTER_PRIMES):
        for c, p in enumerate(MONSTER_PRIMES):
            if p == l:
                continue
            prof = q_profile(p, l, maxdeg=2)
            mat[r][c] = prof.get(2, 0)
            lin[r][c] = prof.get(1, 0)
        if verbose:
            print(f'  l = {l:2d} row done ({time.time() - t0:.1f}s)')
    mism = [(MONSTER_PRIMES[c], MONSTER_PRIMES[r], mat[r][c], NB_QMAT[r][c])
            for r in range(n) for c in range(n)
            if r != c and mat[r][c] != NB_QMAT[r][c]]
    asym = [(MONSTER_PRIMES[c], MONSTER_PRIMES[r], mat[r][c], mat[c][r])
            for r in range(n) for c in range(r + 1, n)
            if mat[r][c] != mat[c][r]]
    if verbose:
        print(f'entries differing from notebook: {len(mism)}')
        for t in mism[:20]:
            print(f'  (p={t[0]}, l={t[1]}): got {t[2]}, notebook {t[3]}')
        print(f'asymmetric pairs in recomputation: {len(asym)}')
        for t in asym[:20]:
            print(f'  Q({t[0]},{t[1]})={t[2]} vs Q({t[1]},{t[0]})={t[3]}')
    return {'qmat': mat, 'linmat': lin, 'mismatches': mism, 'asymmetries': asym}


##########################################################################
# Stage 2/3: kernel-level model of the supersingular locus of X_0(l)^+   #
##########################################################################
# Directed edges of the char-p supersingular l-isogeny graph are pairs
# (j, K) with j in F_{p^2} supersingular and K the kernel polynomial (over
# F_{p^2}: Frobenius^2 = [-p] on E[l], so every kernel poly descends) of an
# order-l subgroup of E_j[l], where E_j is the standard model
#     y^2 = x^3 + 3c x + 2c,  c = j/(1728 - j)
# (E_0: y^2 = x^3 + 1, E_1728: y^2 = x^3 + x).  The model is Frobenius-
# equivariant, so the Frobenius action on edges is literal coefficient
# conjugation.  Points of X_0(l) = edges; points of X_0(l)^+ = orbits under
# the Atkin-Lehner flip w (Velu quotient + dual kernel).  Q counts
# Frobenius 2-orbits of X_0(l)^+ points on the ss locus.
#
# F_{p^2} = F_p[i]/(i^2 - s), s the smallest nonresidue.  Scalars are int
# pairs (u, v); polynomials are np.int64 arrays of shape (2, n), rows =
# (u-part, v-part), low-to-high degree.  Zero polynomial has shape (2, 0).

from functools import lru_cache

# ---------- F_q scalar arithmetic ----------

def smallest_nonresidue(p):
    for s in range(2, p):
        if pow(s, (p - 1) // 2, p) == p - 1:
            return s
    raise ValueError(p)

def fq_add(x, y, p):
    return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)

def fq_sub(x, y, p):
    return ((x[0] - y[0]) % p, (x[1] - y[1]) % p)

def fq_mul(x, y, p, s):
    return ((x[0] * y[0] + s * x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)

def fq_inv(x, p, s):
    nrm = (x[0] * x[0] - s * x[1] * x[1]) % p
    ninv = pow(nrm, p - 2, p)
    return ((x[0] * ninv) % p, (-x[1] * ninv) % p)

def fq_conj(x, p):
    return (x[0] % p, (-x[1]) % p)

def fq_pow(x, e, p, s):
    result = (1, 0)
    while e:
        if e & 1:
            result = fq_mul(result, x, p, s)
        x = fq_mul(x, x, p, s)
        e >>= 1
    return result

# ---------- F_q[x] toolkit ----------

def qzero():
    return np.zeros((2, 0), dtype=np.int64)

def qconst(c, p):
    return qnorm(np.array([[c[0]], [c[1]]], dtype=np.int64), p)

def qnorm(a, p):
    a = np.asarray(a, dtype=np.int64) % p
    nz = np.nonzero(a.any(axis=0))[0]
    return a[:, : nz[-1] + 1] if nz.size else a[:, :0]

def qdeg(a):
    return a.shape[1] - 1

def qlc(a):
    return (int(a[0, -1]), int(a[1, -1]))

def qadd(a, b, p):
    n = max(a.shape[1], b.shape[1])
    out = np.zeros((2, n), dtype=np.int64)
    out[:, : a.shape[1]] += a
    out[:, : b.shape[1]] += b
    return qnorm(out, p)

def qsub(a, b, p):
    n = max(a.shape[1], b.shape[1])
    out = np.zeros((2, n), dtype=np.int64)
    out[:, : a.shape[1]] += a
    out[:, : b.shape[1]] -= b
    return qnorm(out, p)

def qmul(a, b, p, s):
    if a.shape[1] == 0 or b.shape[1] == 0:
        return qzero()
    r = (np.convolve(a[0], b[0]) + s * np.convolve(a[1], b[1])) % p
    im = (np.convolve(a[0], b[1]) + np.convolve(a[1], b[0])) % p
    return qnorm(np.stack([r, im]), p)

def qscale(a, c, p, s):
    return qnorm(np.stack([(a[0] * c[0] + s * a[1] * c[1]),
                           (a[0] * c[1] + a[1] * c[0])]), p)

def qmonic(a, p, s):
    if a.shape[1] == 0:
        return a
    return qscale(a, fq_inv(qlc(a), p, s), p, s)

def qdivmod(a, b, p, s):
    if b.shape[1] == 0:
        raise ZeroDivisionError
    da, db = qdeg(a), qdeg(b)
    if da < db:
        return qzero(), a
    a = a.copy()
    binv = fq_inv(qlc(b), p, s)
    quo = np.zeros((2, da - db + 1), dtype=np.int64)
    for i in range(da - db, -1, -1):
        c = fq_mul((int(a[0, i + db]), int(a[1, i + db])), binv, p, s)
        if c != (0, 0):
            quo[:, i] = c
            seg = a[:, i : i + db + 1]
            seg[0] = (seg[0] - c[0] * b[0] - s * c[1] * b[1]) % p
            seg[1] = (seg[1] - c[0] * b[1] - c[1] * b[0]) % p
    return qnorm(quo, p), qnorm(a, p)

def qmod(a, b, p, s):
    return qdivmod(a, b, p, s)[1]

def qgcd(a, b, p, s):
    while b.shape[1]:
        a, b = b, qmod(a, b, p, s)
    return qmonic(a, p, s)

def qinvmod(v, h, p, s):
    """Inverse of v modulo h (h irreducible or gcd(v,h)=1), extended Euclid."""
    r0, r1 = h, qmod(v, h, p, s)
    s0, s1 = qzero(), qconst((1, 0), p)
    while r1.shape[1]:
        quo, rem = qdivmod(r0, r1, p, s)
        r0, r1 = r1, rem
        s0, s1 = s1, qsub(s0, qmul(quo, s1, p, s), p)
    if qdeg(r0) != 0:
        raise ArithmeticError('qinvmod: not invertible')
    return qmod(qscale(s0, fq_inv(qlc(r0), p, s), p, s), h, p, s)

def qmulmod(a, b, f, p, s):
    return qmod(qmul(a, b, p, s), f, p, s)

def qpowmod(a, e, f, p, s):
    result = qconst((1, 0), p)
    a = qmod(a, f, p, s)
    while e:
        if e & 1:
            result = qmulmod(result, a, f, p, s)
        a = qmulmod(a, a, f, p, s)
        e >>= 1
    return result

def qderiv(a, p):
    if a.shape[1] <= 1:
        return qzero()
    return qnorm(a[:, 1:] * np.arange(1, a.shape[1], dtype=np.int64), p)

def qX(p):
    return qnorm(np.array([[0, 1], [0, 0]], dtype=np.int64), p)

def qeval(a, c, p, s):
    """Evaluate a at the scalar c in F_q (Horner)."""
    acc = (0, 0)
    for i in range(a.shape[1] - 1, -1, -1):
        acc = fq_mul(acc, c, p, s)
        acc = fq_add(acc, (int(a[0, i]), int(a[1, i])), p)
    return acc

# ---------- Brent-Kung modular composition ----------

def qcompose(u, X, f, p, s):
    """u(X) mod f for u, X reduced mod f.  Baby-step giant-step."""
    n = qdeg(f)
    if u.shape[1] == 0:
        return qzero()
    m = int(n ** 0.5) + 1
    baby0 = np.zeros((m, n), dtype=np.int64)
    baby1 = np.zeros((m, n), dtype=np.int64)
    cur = qconst((1, 0), p)
    for jj in range(m):
        baby0[jj, : cur.shape[1]] = cur[0]
        baby1[jj, : cur.shape[1]] = cur[1]
        if jj < m - 1:
            cur = qmulmod(cur, X, f, p, s)
    giant = qmulmod(cur, X, f, p, s)          # X^m mod f
    nch = (u.shape[1] + m - 1) // m
    C0 = np.zeros((nch, m), dtype=np.int64)
    C1 = np.zeros((nch, m), dtype=np.int64)
    for i in range(nch):
        chunk = u[:, i * m : (i + 1) * m]
        C0[i, : chunk.shape[1]] = chunk[0]
        C1[i, : chunk.shape[1]] = chunk[1]
    V0 = (C0 @ baby0 + s * (C1 @ baby1)) % p
    V1 = (C0 @ baby1 + C1 @ baby0) % p
    result = qzero()
    for i in range(nch - 1, -1, -1):
        result = qmulmod(result, giant, f, p, s)
        result = qadd(result, qnorm(np.stack([V0[i], V1[i]]), p), p)
    return result

# ---------- equal-degree factor extraction ----------

def _q_random_poly(deg, p, rng):
    return qnorm(rng.integers(0, p, size=(2, deg)), p)

def one_irreducible_factor(f, r, Xq, p, s, rng):
    """One irreducible factor of squarefree f; all factors have degree r."""
    n = qdeg(f)
    assert n % r == 0, (n, r)
    if n == r:
        return qmonic(f, p, s)
    q = p * p
    for _ in range(60):
        u = _q_random_poly(n, p, rng)
        # trace to the F_q-subalgebra: t = sum u^(q^i), i < r
        t = u
        cur = u
        for _ in range(r - 1):
            cur = qcompose(cur, Xq, f, p, s)
            t = qadd(t, cur, p)
        for probe in (qpowmod(t, (q - 1) // 2, f, p, s), t):
            for shift in ((1, 0), (p - 1, 0), (0, 0)):
                g = qgcd(qsub(probe, qconst(shift, p), p) if shift != (0, 0)
                         else probe, f, p, s)
                dg = qdeg(g)
                if 0 < dg < n:
                    piece = g if dg <= n - dg else qdivmod(f, g, p, s)[0]
                    return one_irreducible_factor(
                        piece, r, qmod(Xq, piece, p, s), p, s, rng)
    raise RuntimeError('one_irreducible_factor: no split found')

# ---------- division polynomials (x-only) ----------

def divpolys_mod(a, b, h, upto, p, s):
    """f_n mod h for n <= upto, where psi_n = f_n * (2y)^(n even) on
    y^2 = x^3 + a x + b.  All arithmetic in F_q[x]/(h)."""
    X = qmod(qX(p), h, p, s)
    aa, bb = qconst(a, p), qconst(b, p)
    W = qmod(qnorm(np.stack([
        np.array([4 * b[0] % p, 4 * a[0] % p, 0, 4]),
        np.array([4 * b[1] % p, 4 * a[1] % p, 0, 0])]), p), h, p, s)
    mm = lambda u, v: qmulmod(u, v, h, p, s)
    f = {0: qzero(), 1: qconst((1, 0), p), 2: qconst((1, 0), p)}
    # f_3 = 3x^4 + 6a x^2 + 12b x - a^2 ; f_4 = 2(x^6 + 5a x^4 + 20b x^3
    #        - 5a^2 x^2 - 4ab x - 8b^2 - a^3)
    x2 = mm(X, X); x3 = mm(x2, X); x4 = mm(x3, X)
    a2 = fq_mul(a, a, p, s); ab = fq_mul(a, b, p, s); b2 = fq_mul(b, b, p, s)
    a3 = fq_mul(a2, a, p, s)
    f[3] = qadd(qadd(qscale(x4, (3, 0), p, s), qscale(x2, fq_mul((6, 0), a, p, s), p, s), p),
                qadd(qscale(X, fq_mul((12, 0), b, p, s), p, s),
                     qconst(fq_mul((p - 1, 0), a2, p, s), p), p), p)
    x6 = mm(x3, x3)
    t4 = qadd(x6, qscale(x4, fq_mul((5, 0), a, p, s), p, s), p)
    t4 = qadd(t4, qscale(x3, fq_mul((20, 0), b, p, s), p, s), p)
    t4 = qsub(t4, qscale(x2, fq_mul((5, 0), a2, p, s), p, s), p)
    t4 = qsub(t4, qscale(X, fq_mul((4, 0), ab, p, s), p, s), p)
    t4 = qsub(t4, qconst(fq_add(fq_mul((8, 0), b2, p, s), a3, p), p), p)
    f[4] = qscale(t4, (2, 0), p, s)
    W2 = mm(W, W)
    for n in range(5, upto + 1):
        if n in f:
            continue
        m = n // 2
        if n % 2 == 1:
            t1 = mm(f[m + 2], mm(mm(f[m], f[m]), f[m]))
            t2 = mm(f[m - 1], mm(mm(f[m + 1], f[m + 1]), f[m + 1]))
            f[n] = qsub(mm(W2, t1), t2, p) if m % 2 == 0 else qsub(t1, mm(W2, t2), p)
        else:
            t1 = mm(f[m + 2], mm(f[m - 1], f[m - 1]))
            t2 = mm(f[m - 2], mm(f[m + 1], f[m + 1]))
            f[n] = mm(f[m], qsub(t1, t2, p))
    return f, W

def psi_l_full(a, b, l, p, s):
    """The full division polynomial f_l (x-only part), NOT reduced."""
    # reuse divpolys_mod with a modulus of degree > deg f_l (i.e. no reduction)
    big = np.zeros((2, (l * l - 1) // 2 + 2), dtype=np.int64)
    big[0, -1] = 1
    f, _ = divpolys_mod(a, b, qnorm(big, p), l, p, s)
    return f[l]

# ---------- linear algebra mod p, minimal polynomials over F_q ----------

def _rowreduce_solve(A, b, p):
    """Solve A x = b mod p (A: 2d int64 array).  Returns x or None."""
    A = A.astype(np.int64) % p
    b = b.astype(np.int64) % p
    nr, nc = A.shape
    M = np.concatenate([A, b.reshape(-1, 1)], axis=1)
    piv_cols = []
    row = 0
    for col in range(nc):
        sel = None
        for rr in range(row, nr):
            if M[rr, col] % p:
                sel = rr
                break
        if sel is None:
            continue
        M[[row, sel]] = M[[sel, row]]
        inv = pow(int(M[row, col]), p - 2, p)
        M[row] = (M[row] * inv) % p
        mask = np.arange(nr) != row
        M[mask] = (M[mask] - np.outer(M[mask, col], M[row])) % p
        piv_cols.append(col)
        row += 1
        if row == nr:
            break
    # consistency
    for rr in range(row, nr):
        if M[rr, -1] % p:
            return None
    x = np.zeros(nc, dtype=np.int64)
    for k, col in enumerate(piv_cols):
        x[col] = M[k, -1]
    return x

def _rh_vec(v, r):
    out = np.zeros(2 * r, dtype=np.int64)
    out[: v.shape[1]] = v[0]
    out[r : r + v.shape[1]] = v[1]
    return out

def minpoly_fq(beta, h, p, s):
    """Minimal polynomial over F_q of beta in F_q[x]/(h) (monic q-poly)."""
    r = qdeg(h)
    pows = [qconst((1, 0), p)]
    for _ in range(r):
        pows.append(qmulmod(pows[-1], beta, h, p, s))
    for m in range(1, r + 1):
        # columns: vec(beta^k), vec(i*beta^k) for k < m
        cols = []
        for k in range(m):
            v = pows[k]
            cols.append(_rh_vec(v, r))
            iv = qnorm(np.stack([(s * v[1]) % p, v[0]]), p)   # i * v
            cols.append(_rh_vec(iv, r))
        A = np.stack(cols, axis=1)
        b = (-_rh_vec(pows[m], r)) % p
        x = _rowreduce_solve(A, b, p)
        if x is not None:
            coefs = np.zeros((2, m + 1), dtype=np.int64)
            for k in range(m):
                coefs[0, k] = x[2 * k]
                coefs[1, k] = x[2 * k + 1]
            coefs[0, m] = 1
            return qnorm(coefs, p)
    raise RuntimeError('minpoly_fq failed')

# ---------- kernel enumeration at a vertex ----------

def _orbit_reps_mod_pm(l, p):
    """Representatives (and orbits) of <p, -1> acting on (Z/l)^* by mult,
    intersected with 1..(l-1)/2."""
    seen = set()
    reps = []
    for k in range(1, (l + 1) // 2):
        if k in seen:
            continue
        reps.append(k)
        orb = set()
        cur = k
        while cur not in orb:
            orb.add(cur)
            orb.add((l - cur) % l)
            cur = (cur * p) % l
        seen |= {min(t, l - t) for t in orb}
    return reps

def frob_order_mod_pm(l, p):
    """Order of p in (Z/l)^* / {+-1} = common degree of psi_l factors."""
    cur = p % l
    for r in range(1, l):
        if cur == 1 or cur == l - 1:
            return r
        cur = (cur * p) % l
    raise ValueError((l, p))

def kernel_from_factor(a, b, h, l, p, s):
    """Kernel polynomial (monic, degree (l-1)/2) of the order-l subgroup
    through a root of the irreducible factor h of psi_l."""
    d = (l - 1) // 2
    f, W = divpolys_mod(a, b, h, (l + 3) // 2, p, s)
    X = qmod(qX(p), h, p, s)
    K = h
    covered = qdeg(h)
    for k in _orbit_reps_mod_pm(l, p):
        if covered >= d:
            break
        if k == 1:
            continue
        num = qmulmod(f[k - 1], f[k + 1], h, p, s)
        den = qmulmod(f[k], f[k], h, p, s)
        if k % 2 == 0:
            den = qmulmod(W, den, h, p, s)
        else:
            num = qmulmod(W, num, h, p, s)
        xi = qsub(X, qmulmod(num, qinvmod(den, h, p, s), h, p, s), p)
        g = minpoly_fq(xi, h, p, s)
        rem = qmod(K, g, p, s)
        if rem.shape[1] != 0:      # new factor of this kernel
            K = qmul(K, g, p, s)
            covered += qdeg(g)
    assert qdeg(K) == d, (qdeg(K), d)
    return qmonic(K, p, s)

def kernels_at_vertex(a, b, l, p, s, seed=0):
    """All l+1 kernel polynomials of E: y^2 = x^3 + ax + b (supersingular)."""
    rng = np.random.default_rng(seed)
    if l == 2:
        cubic = qnorm(np.stack([np.array([b[0], a[0], 0, 1]),
                                np.array([b[1], a[1], 0, 0])]), p)
        roots = qroots(cubic, p, s)
        assert len(roots) == 3, roots
        return [qnorm(np.stack([np.array([(-x0[0]) % p, 1]),
                                np.array([(-x0[1]) % p, 0])]), p) for x0 in roots]
    psi = qmonic(psi_l_full(a, b, l, p, s), p, s)
    assert qdeg(psi) == (l * l - 1) // 2
    r = frob_order_mod_pm(l, p)
    Xq = qpowmod(qX(p), p * p, psi, p, s)
    kernels = []
    remaining = psi
    while qdeg(remaining) > 0:
        h = one_irreducible_factor(remaining, r, qmod(Xq, remaining, p, s), p, s, rng)
        K = kernel_from_factor(a, b, h, l, p, s)
        quo, rem = qdivmod(remaining, K, p, s)
        assert rem.shape[1] == 0, 'kernel does not divide psi_l'
        remaining = quo
        kernels.append(K)
    assert len(kernels) == l + 1, (len(kernels), l + 1)
    return kernels

def qroots(f, p, s):
    """All roots in F_q of f (brute-force vectorized scan; p smallish)."""
    U, V = _fq_grid(p)
    acc0 = np.zeros_like(U)
    acc1 = np.zeros_like(U)
    for i in range(f.shape[1] - 1, -1, -1):
        acc0, acc1 = (acc0 * U + s * acc1 * V) % p, (acc0 * V + acc1 * U) % p
        acc0 = (acc0 + int(f[0, i])) % p
        acc1 = (acc1 + int(f[1, i])) % p
    hit = (acc0 == 0) & (acc1 == 0)
    return [(int(u), int(v)) for u, v in zip(U[hit], V[hit])]

@lru_cache(maxsize=8)
def _fq_grid(p):
    U, V = np.meshgrid(np.arange(p, dtype=np.int64), np.arange(p, dtype=np.int64))
    return U.ravel(), V.ravel()

# ---------- standard models, vertices, modular-polynomial checks ----------

def std_model(j, p, s):
    """Frobenius-equivariant model of j: y^2 = x^3 + 3c x + 2c,
    c = j/(1728-j); (0,1) at j=0 and (1,0) at j=1728."""
    if j == (0, 0):
        return (0, 0), (1, 0)
    if j == (1728 % p, 0):
        return (1, 0), (0, 0)
    c = fq_mul(j, fq_inv(fq_sub((1728 % p, 0), j, p), p, s), p, s)
    return fq_mul((3, 0), c, p, s), fq_mul((2, 0), c, p, s)

_phi_cache = {}

def _load_phi(l):
    if l not in _phi_cache:
        with open(_DATA_DIR / 'classical_modpolys.json') as f:
            raw = json.load(f)
        if str(l) not in raw:
            raise ValueError(f'no classical modular polynomial for l = {l}')
        _phi_cache[l] = raw[str(l)]
    return _phi_cache[l]

def phi_check(l, j1, j2, p, s):
    """Phi_l(j1, j2) == 0 over F_q."""
    M = _load_phi(l)
    n = len(M)
    pows1 = [(1, 0)]
    pows2 = [(1, 0)]
    for _ in range(n - 1):
        pows1.append(fq_mul(pows1[-1], j1, p, s))
        pows2.append(fq_mul(pows2[-1], j2, p, s))
    acc = (0, 0)
    for i in range(n):
        for k in range(n):
            m = M[i][k] % p
            if m:
                acc = fq_add(acc, fq_mul((m, 0), fq_mul(pows1[i], pows2[k], p, s), p, s), p)
    return acc == (0, 0)

def phi_fiber_poly(l, j, p, s):
    """Phi_l(j, Y) as a q-poly in Y."""
    M = _load_phi(l)
    n = len(M)
    pows = [(1, 0)]
    for _ in range(n - 1):
        pows.append(fq_mul(pows[-1], j, p, s))
    coefs = np.zeros((2, n), dtype=np.int64)
    for k in range(n):
        acc = (0, 0)
        for i in range(n):
            m = M[i][k] % p
            if m:
                acc = fq_add(acc, fq_mul((m, 0), pows[i], p, s), p)
        coefs[:, k] = acc
    return qnorm(coefs, p)

def ss_vertex_set(p, s):
    """All supersingular j in F_{p^2}: 2-isogeny closure of the rational ones."""
    seeds = [(j % p, 0) for j in rational_ss_js(p)]
    seen = set(seeds)
    stack = list(seeds)
    while stack:
        j = stack.pop()
        for j2 in qroots(phi_fiber_poly(2, j, p, s), p, s):
            if j2 not in seen:
                seen.add(j2)
                stack.append(j2)
    expected = p // 12 + {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]
    assert len(seen) == expected, (p, len(seen), expected)
    return sorted(seen)

# ---------- automorphism normalization of kernels ----------

@lru_cache(maxsize=64)
def _cube_roots_of_unity(p, s):
    """Primitive cube roots of unity in F_q (p != 3)."""
    quad = qnorm(np.array([[1, 1, 1], [0, 0, 0]], dtype=np.int64), p)
    return tuple(qroots(quad, p, s))

def kernel_tuple(K):
    return tuple((int(K[0, i]), int(K[1, i])) for i in range(K.shape[1]))

def kernel_canonical(K, j, p, s):
    """Canonical form of the kernel under Aut(E_j): x -> zeta*x scalings."""
    d = qdeg(K)
    scalars = [(1, 0)]
    if j == (0, 0):
        scalars += list(_cube_roots_of_unity(p, s))
    elif j == (1728 % p, 0):
        scalars.append((p - 1, 0))
    best = None
    for z in scalars:
        zp = (1, 0)
        Kz = K.copy()
        for i in range(d - 1, -1, -1):        # coefficient of x^i scales by z^(d-i)
            zp = fq_mul(zp, z, p, s)
            Kz[:, i] = fq_mul((int(K[0, i]), int(K[1, i])), zp, p, s)
        t = kernel_tuple(qnorm(Kz, p))
        if best is None or t < best:
            best = t
    return best

# ---------- edges and the two involutions ----------

def kernels_with_seeds(a, b, l, p, s, seed=0):
    """[(kernel poly, irreducible seed factor)] for all l+1 subgroups."""
    rng = np.random.default_rng(seed)
    if l == 2:
        ks = kernels_at_vertex(a, b, l, p, s, seed)
        return [(K, K) for K in ks]
    psi = qmonic(psi_l_full(a, b, l, p, s), p, s)
    assert qdeg(psi) == (l * l - 1) // 2
    r = frob_order_mod_pm(l, p)
    Xq = qpowmod(qX(p), p * p, psi, p, s)
    out = []
    remaining = psi
    while qdeg(remaining) > 0:
        h = one_irreducible_factor(remaining, r, qmod(Xq, remaining, p, s), p, s, rng)
        K = kernel_from_factor(a, b, h, l, p, s)
        quo, rem = qdivmod(remaining, K, p, s)
        assert rem.shape[1] == 0
        remaining = quo
        out.append((K, h))
    assert len(out) == l + 1
    return out

def _powersums(K, p, s):
    """p1, p2, p3 of the roots of monic K (elementary symmetric via coeffs)."""
    d = qdeg(K)
    c = lambda i: (int(K[0, i]), int(K[1, i])) if 0 <= i < K.shape[1] else (0, 0)
    e1 = fq_sub((0, 0), c(d - 1), p)
    e2 = c(d - 2) if d >= 2 else (0, 0)
    e3 = fq_sub((0, 0), c(d - 3), p) if d >= 3 else (0, 0)
    p1 = e1
    p2 = fq_sub(fq_mul(e1, p1, p, s), fq_mul((2, 0), e2, p, s), p)
    p3 = fq_add(fq_sub(fq_mul(e1, p2, p, s), fq_mul(e2, p1, p, s), p),
                fq_mul((3, 0), e3, p, s), p)
    return p1, p2, p3

def velu_codomain(a, b, K, l, p, s):
    """(a', b') of E/C for kernel poly K (Velu, full-set sums)."""
    d = qdeg(K)
    tau = (1, 0) if l == 2 else (2, 0)
    p1, p2, p3 = _powersums(K, p, s)
    dd = (d % p, 0)
    # Velu with S a half-set: t_Q = tau*(3x^2+a), u_Q = 4(x^3+ax+b)
    # t = sum t_Q ; w = sum u_Q + x_Q t_Q  (u_Q = 0 identically for l = 2)
    t = fq_mul(tau, fq_add(fq_mul((3, 0), p2, p, s), fq_mul(dd, a, p, s), p), p, s)
    w = fq_mul(tau, fq_add(fq_mul((3, 0), p3, p, s), fq_mul(a, p1, p, s), p), p, s)
    if l != 2:
        u = fq_add(fq_add(p3, fq_mul(a, p1, p, s), p), fq_mul(b, dd, p, s), p)
        w = fq_add(w, fq_mul((4, 0), u, p, s), p)
    av = fq_sub(a, fq_mul((5, 0), t, p, s), p)
    bv = fq_sub(b, fq_mul((7, 0), w, p, s), p)
    return av, bv

def fq_j_invariant(a, b, p, s):
    a3 = fq_mul((4, 0), fq_mul(fq_mul(a, a, p, s), a, p, s), p, s)
    den = fq_add(a3, fq_mul((27, 0), fq_mul(b, b, p, s), p, s), p)
    return fq_mul(fq_mul((1728 % p, 0), a3, p, s), fq_inv(den, p, s), p, s)

def phi_x_at(a, b, K, l, h, p, s):
    """phi_x(alpha) in R_h = F_q[x]/(h), alpha = class of x, for the Velu
    isogeny with kernel K (alpha NOT a root of K)."""
    d = qdeg(K)
    tau = 1 if l == 2 else 2
    p1, _, _ = _powersums(K, p, s)
    X = qmod(qX(p), h, p, s)
    Kh = qmod(K, h, p, s)
    K1 = qmod(qderiv(K, p), h, p, s)
    K2 = qmod(qderiv(qderiv(K, p), p), h, p, s)
    inv = qinvmod(Kh, h, p, s)
    S1 = qmulmod(K1, inv, h, p, s)
    S2 = qmulmod(qsub(qmulmod(K1, K1, h, p, s), qmulmod(Kh, K2, h, p, s), p),
                 qmulmod(inv, inv, h, p, s), h, p, s)
    X2 = qmulmod(X, X, h, p, s)
    X3 = qmulmod(X2, X, h, p, s)
    aC, bC = qconst(a, p), qconst(b, p)
    W0 = qadd(qadd(X3, qmulmod(aC, X, h, p, s), p), bC, p)        # x^3+ax+b
    W0p = qadd(qscale(X2, (3, 0), p, s), aC, p)                    # 3x^2+a
    dX = qscale(X, (d % p, 0), p, s)
    p1C = qconst(p1, p)
    # sum t_Q/(alpha-x) = 3(alpha^2 S1 - d*alpha - p1) + a S1
    tsum = qadd(qscale(qsub(qsub(qmulmod(X2, S1, h, p, s), dX, p), p1C, p), (3, 0), p, s),
                qmulmod(aC, S1, h, p, s), p)
    # sum u_Q/(alpha-x)^2 = 4[W0(a)S2 - W0'(a)S1 + 3 d alpha - (d alpha - p1)]
    usum = qsub(qmulmod(W0, S2, h, p, s), qmulmod(W0p, S1, h, p, s), p)
    usum = qadd(usum, qscale(dX, (3, 0), p, s), p)
    usum = qsub(usum, qsub(dX, p1C, p), p)
    usum = qscale(usum, (4, 0), p, s)
    # t_Q carries the half-set doubling (tau); u_Q enters once per Q in S
    return qadd(X, qadd(qscale(tsum, (tau % p, 0), p, s), usum, p), p)

def dual_kernel_match(Khat, av, bv, jt, kernels_t, p, s):
    """Index of the kernel at target vertex jt matching Khat (kernel of the
    dual, in Velu-codomain coordinates).  Matching is up to the scaling
    x -> u^2 x with u^4 = av/a', u^6 = bv/b'; at jt in {0, 1728} the leftover
    root ambiguity is exactly Aut(E_jt), which the edge quotient kills."""
    at, bt = std_model(jt, p, s)
    d = qdeg(Khat)
    ch = lambda i: (int(Khat[0, i]), int(Khat[1, i])) if i < Khat.shape[1] else (0, 0)
    matches = []
    for idx, Kt in enumerate(kernels_t):
        if qdeg(Kt) != d:
            continue
        ct = lambda i: (int(Kt[0, i]), int(Kt[1, i])) if i < Kt.shape[1] else (0, 0)
        if jt == (0, 0):
            u6 = bv                                   # b' = 1
            ok = all(
                (ch(i) == (0, 0)) == (ct(i) == (0, 0)) and
                (ch(i) == (0, 0) or
                 fq_pow(ch(i), 3, p, s) ==
                 fq_mul(fq_pow(u6, d - i, p, s), fq_pow(ct(i), 3, p, s), p, s))
                for i in range(d))
        elif jt == (1728 % p, 0):
            u4 = av                                   # a' = 1
            ok = all(
                (ch(i) == (0, 0)) == (ct(i) == (0, 0)) and
                (ch(i) == (0, 0) or
                 fq_pow(ch(i), 2, p, s) ==
                 fq_mul(fq_pow(u4, d - i, p, s), fq_pow(ct(i), 2, p, s), p, s))
                for i in range(d))
        else:
            u2 = fq_mul(fq_mul(bv, at, p, s), fq_inv(fq_mul(av, bt, p, s), p, s), p, s)
            ok = all(ch(i) == fq_mul(fq_pow(u2, d - i, p, s), ct(i), p, s)
                     for i in range(d))
        if ok:
            matches.append(idx)
    return matches

# ---------- graph assembly and Q ----------

def build_edge_model(p, l, seed=0, verbose=False):
    """All edges (points of X_0(l) on the ss locus over F_{p^2}) with the
    Frobenius and Atkin-Lehner involutions as permutations."""
    if p == l:
        raise ValueError('p and l must be distinct')
    s = smallest_nonresidue(p)
    vertices = ss_vertex_set(p, s)
    vset = set(vertices)
    edges = []                       # [{j, K, seed_h}]
    index = {}                       # (j, canonical) -> edge id
    kern_lists = {}                  # j -> [(K, h, edge_id)]
    t0 = time.time()
    for j in vertices:
        a, b = std_model(j, p, s)
        ks = kernels_with_seeds(a, b, l, p, s, seed=seed)
        kern_lists[j] = []
        for K, h in ks:
            key = (j, kernel_canonical(K, j, p, s))
            if key not in index:
                index[key] = len(edges)
                edges.append({'j': j, 'K': K, 'h': h})
            kern_lists[j].append((K, h, index[key]))
        if verbose:
            print(f'    vertex {j}: {len(ks)} kernels ({time.time()-t0:.1f}s)')
    # Frobenius
    F = []
    for e in edges:
        jt = fq_conj(e['j'], p)
        Kc = qnorm(np.stack([e['K'][0], (-e['K'][1]) % p]), p)
        F.append(index[(jt, kernel_canonical(Kc, jt, p, s))])
    # Atkin-Lehner
    Wp = []
    for eid, e in enumerate(edges):
        j, K = e['j'], e['K']
        a, b = std_model(j, p, s)
        av, bv = velu_codomain(a, b, K, l, p, s)
        jt = fq_j_invariant(av, bv, p, s)
        assert jt in vset, (j, jt)
        assert phi_check(l, j, jt, p, s), (p, l, j, jt)
        # a root of any other kernel's seed gives beta = phi_x(alpha) in C-dual
        # (any subgroup other than C itself works, even within the same
        # Aut-orbit/edge — alpha just must not be a root of K)
        Kt = kernel_tuple(K)
        h = next(hh for KK, hh, ii in kern_lists[j] if kernel_tuple(KK) != Kt)
        beta = phi_x_at(a, b, K, l, h, p, s)
        if l == 2:
            Khat = qnorm(np.stack([np.array([(-beta[0, 0]) % p if beta.shape[1] else 0, 1]),
                                   np.array([(-beta[1, 0]) % p if beta.shape[1] else 0, 0])]), p)
        else:
            hpp = minpoly_fq(beta, h, p, s)
            Khat = kernel_from_factor(av, bv, hpp, l, p, s)
        kt = [KK for KK, hh, ii in kern_lists[jt]]
        ids = [ii for KK, hh, ii in kern_lists[jt]]
        midx = dual_kernel_match(Khat, av, bv, jt, kt, p, s)
        hit_ids = {ids[i] for i in midx}
        assert len(hit_ids) == 1, f'dual match not unique: edges {hit_ids} at {jt}'
        Wp.append(hit_ids.pop())
    # structural checks
    n = len(edges)
    assert all(F[F[i]] == i for i in range(n)), 'F not an involution'
    assert all(Wp[Wp[i]] == i for i in range(n)), 'w not an involution'
    assert all(F[Wp[i]] == Wp[F[i]] for i in range(n)), 'F and w do not commute'
    return {'p': p, 'l': l, 's': s, 'edges': edges, 'F': F, 'w': Wp}

def count_Q(model):
    """Q_full and Q_slice: Frobenius 2-orbits of X_0(l)^+ points on the ss
    locus (all, resp. those whose j-pair meets the F_p-rational vertices)."""
    p, F, w, edges = model['p'], model['F'], model['w'], model['edges']
    n = len(edges)
    seen = [False] * n
    quad_full = quad_slice = pts_full = pts_rat = 0
    for i in range(n):
        if seen[i]:
            continue
        orb = {i, w[i]}
        for t in orb:
            seen[t] = True
        pts_full += 1
        quadratic = F[i] not in orb
        touches_rational = any(edges[t]['j'][1] == 0 for t in orb)
        if quadratic:
            quad_full += 1
            if touches_rational:
                quad_slice += 1
        elif touches_rational:
            pts_rat += 1
    assert quad_full % 2 == 0 and quad_slice % 2 == 0
    return {'Q_full': quad_full // 2, 'Q_slice': quad_slice // 2,
            'points': pts_full, 'n_edges': len(edges)}

def stage2_Q(p, l, verbose=False):
    return count_Q(build_edge_model(p, l, verbose=verbose))

def involution_profile(p, l):
    """Fixed-point census of the Klein group <F, w> on the ss points of
    X_0(l) in char p, and the resulting closed form for Q_full.

    D = # ss points of X_0(l) (moduli, Aut-orbits merged)
    A = # Frobenius-fixed points   (F_p-rational isogenies, incl. isogenies
        to the quadratic twist — 'potential endomorphisms')
    B = # points with F(e) = w(e)  ('Hermitian': conjugate-dual = itself)
    C = # w-fixed points           (self-dual isogenies, CM by sqrt(-l))

    Then  #points(X_0(l)^+, ss) = (D + C)/2,  #rational = (A + B)/2,
    and   Q_full = (D + C - A - B)/4."""
    m = build_edge_model(p, l)
    F, w = m['F'], m['w']
    n = len(m['edges'])
    D = n
    A = sum(F[i] == i for i in range(n))
    B = sum(F[i] == w[i] for i in range(n))
    C = sum(w[i] == i for i in range(n))
    q = count_Q(m)
    assert (D + C - A - B) % 4 == 0
    assert (D + C - A - B) // 4 == q['Q_full'], (D, A, B, C, q)
    return {'p': p, 'l': l, 'D': D, 'A': A, 'B': B, 'C': C,
            'points': (D + C) // 2, 'rational': (A + B) // 2,
            'Q_full': q['Q_full'], 'Q_slice': q['Q_slice']}

def stage2_calibrate(pairs=None, verbose=True):
    """Check Q_slice against the stage-1 (notebook) computation."""
    if pairs is None:
        small = [5, 7, 11, 13]
        pairs = [(p, l) for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]
                 for l in small if p != l]
    bad = []
    for p, l in pairs:
        t0 = time.time()
        got = stage2_Q(p, l)
        want = Q(p, l)
        tag = 'OK' if got['Q_slice'] == want else 'MISMATCH'
        if tag == 'MISMATCH':
            bad.append((p, l, got['Q_slice'], want))
        if verbose:
            print(f"  (p={p:3d}, l={l:2d}): stage2 Q_slice={got['Q_slice']:3d} "
                  f"stage1 Q={want:3d}  Q_full={got['Q_full']:3d}  "
                  f"[{got['n_edges']} edges, {time.time()-t0:.1f}s] {tag}")
    return bad


#########################################################
# Genus of X_0(N)* and the conjecture Q = g* - g(l)+    #
#########################################################

from identities import clgr_size_gen as _h

def _fix_wQ(Q, others, restrict=True):
    """# fixed points of w_Q on X_0(Q * prod(others)), Q >= 5 prime."""
    loc4 = 1
    locQ = 1
    for q in others:
        loc4 *= 1 + quad_rec(-4 * Q, q)
        locQ *= 1 + quad_rec(-Q, q)
    n = _h(-4 * Q) * loc4
    if (-Q) % 4 == 1:
        n += _h(-Q) * locQ
    return n

def genus_x0(N_primes):
    """Genus of X_0(N), N = squarefree product of the given primes."""
    mu = 1
    nu2 = 1
    nu3 = 1
    for q in N_primes:
        mu *= q + 1
        nu2 *= 1 + quad_rec(-4, q)
        nu3 *= 1 + quad_rec(-3, q)
    nuinf = 2 ** len(N_primes)
    g12 = 12 + mu - 3 * nu2 - 4 * nu3 - 6 * nuinf
    assert g12 % 12 == 0
    return g12 // 12

def genus_x0_star(p, l):
    """Genus of X_0(pl)* = X_0(pl)/<w_p, w_l>, p, l >= 5 distinct primes."""
    g0 = genus_x0((p, l))
    nu = _fix_wQ(p, (l,)) + _fix_wQ(l, (p,)) + _fix_wQ(p * l, ())
    gs = 1 + (g0 - 1) / 4 - nu / 8
    assert gs == int(gs) and gs >= 0, (p, l, g0, nu, gs)
    return int(gs)

def genus_x0_plus(l):
    """Genus of X_0(l)+, l >= 5 prime."""
    g0 = genus_x0((l,))
    nu = _fix_wQ(l, ())
    gp = 1 + (g0 - 1) / 2 - nu / 4
    assert gp == int(gp) and gp >= 0, (l, g0, nu, gp)
    return int(gp)

def _fix_wQ_needs(Q):
    return Q >= 5

def test_genus_conjecture(verbose=True):
    """Test Q_full(p,l) == g(X_0(pl)*) - g(X_0(l)+) against everything:
    the notebook matrix (Monster p: Q_full = notebook Q) and the stage-4
    cache.  Restricted to p, l >= 5 (w_2, w_3 fixed-point conventions)."""
    bad = []
    tested = 0
    # notebook matrix: characteristic Monster => Q_full = Q_notebook
    M = MONSTER_PRIMES
    for r, l in enumerate(M):
        for c, p in enumerate(M):
            if p == l or p < 5 or l < 5:
                continue
            pred = genus_x0_star(p, l) - genus_x0_plus(l)
            tested += 1
            if pred != NB_QMAT[r][c]:
                bad.append(('NB', p, l, NB_QMAT[r][c], pred))
    # stage-4 cache (includes non-Monster p and l)
    res = _results_load()
    for key, v in res.items():
        p, l = (int(x) for x in key.split(','))
        if p < 5 or l < 5:
            continue
        pred = genus_x0_star(p, l) - genus_x0_plus(l)
        tested += 1
        if pred != v['Q_full']:
            bad.append(('S4', p, l, v['Q_full'], pred))
    if verbose:
        print(f'tested {tested} values; {len(bad)} mismatches')
        for t in bad[:30]:
            print('  ', t)
    return bad

from nt import quad_rec

##########################################################################
# Per-discriminant pairing verification (for the identity paper)        #
##########################################################################
# The identity x0l_fp_card(p,l) + supsingtrace(l,p) = 2p+2 should hold
# discriminant by discriminant: for each trace a in [0, 2 sqrt p], with
# d = a^2 - 4p, the a-contribution to #X_0(l)(F_p) plus the a-contribution
# to Tr(char-l degree-p graph) should equal 2 * clgr_sum(d), up to explicit
# corrections at j = 0, 1728 (the only j's straddling several classes,
# since only O_{-3}, O_{-4} have extra units) and at ramified d.

from identities import (clgr_sum as _clgr_sum, clgr_size_gen as _clgr_gen,
                        supsingtrace as _sstr, x0l_fp_card as _x0card,
                        disc_closure as _disc_closure)
from nt import primesBetween, discfac as _discfac

def identity_sweep(N=200, verbose=True):
    """Check the 2p+2 identity for all ordered pairs of primes 5 <= p,l < N,
    plus the forced anchors: supsingtrace(l, p) == p+1 for l-1 | 12."""
    ps = [q for q in primesBetween(2, N)]
    bad = []
    for p in ps:
        for l in ps:
            if p == l:
                continue
            if _x0card(p, l) + _sstr(l, p) != 2 * p + 2:
                bad.append(('identity', p, l))
        for l in (2, 3, 5, 7, 13):
            if l != p and _sstr(l, p) != p + 1:
                bad.append(('anchor', p, l))
    if verbose:
        n = len(ps) * (len(ps) - 1)
        print(f'identity checked on {n} ordered pairs (p,l < {N}): '
              f'{len(bad)} failures')
        for t in bad[:20]:
            print('  ', t)
    return bad

def _x0_contrib_by_a(p, l):
    """Noncuspidal #X_0(l)(F_p) decomposed by trace a >= 0 (mirrors
    x0l_fp_card / ecfp_disc_l_isocts, scoped to one isogeny class).

    The orders O_{-3}, O_{-4} (extra units) are the only ones straddling
    several classes (j = 0, 1728 quartic/sextic twists); their contribution
    is attached to the smallest a where they appear — the same convention
    as the d3seen/d4seen dedup in identities.supsingtrace."""
    out = {}
    a = 1
    while a * a < 4 * p:
        cl = _disc_closure(a * a - 4 * p)
        clset = set(cl)
        total = 0
        for d0 in cl:
            if d0 in (-3, -4):
                # j = 0, 1728 straddle several classes (twists with distinct
                # traces); their X_0-points aggregate kernels across twists
                # and are verified in the special block, not per class
                continue
            qr = quad_rec(d0, l)
            if d0 * l * l in clset:            # vertices with descending edges
                if d0 == -3:
                    nu = 1 + qr + (l - qr) // 3
                elif d0 == -4:
                    nu = 1 + qr + (l - qr) // 2
                else:
                    nu = l + 1
            else:
                nu = 1 + qr
            total += nu * _clgr_gen(d0)
        out[a] = total
        a += 1
    # supersingular classes (a = 0), as in x0l_fp_card
    out[0] = _clgr_sum(-4 * p) if quad_rec(-p, l) == 1 else 0
    return out

def _trace_contrib_by_a(l, p):
    """supsingtrace(l, p) decomposed by a; the O_{-3}, O_{-4} parts are
    excluded per class (verified in the special block instead)."""
    out = {}
    d3seen = d4seen = False
    a = 0
    while a * a < 4 * p:
        term = 0
        d, c = _discfac(a * a - 4 * p)
        if quad_rec(d, l) < 1:
            while c % l == 0:
                c //= l
            h = _clgr_sum(d * c * c)
            ram = (d % l == 0) or (d % p == 0)
            mult = 1 if ram else 2
            if d == -3:
                term = mult * (h if not d3seen else h - 1)
                d3seen = True
            elif d == -4:
                term = mult * (h if not d4seen else h - 1)
                d4seen = True
            else:
                term = mult * h
        out[a] = term
        a += 1
    return out

def pairing_table(p, l, verbose=True):
    """Per-a table: X0(a), Tr(a), 2*clgr_sum(a^2-4p), delta.  Verifies the
    global identity and catalogs where the per-a pairing needs corrections."""
    x0 = _x0_contrib_by_a(p, l)
    tr = _trace_contrib_by_a(l, p)
    all_ds = set()
    fams = {-3: [], -4: []}
    for a in sorted(set(x0) | set(tr)):
        d = a * a - 4 * p
        all_ds |= set(_disc_closure(d))
        f = _discfac(d)[0]
        if a >= 1 and f in fams:
            fams[f].append(a)
    rows = []
    for a in sorted(set(x0) | set(tr)):
        d = a * a - 4 * p
        f = _discfac(d)[0]
        if a >= 1 and f in fams:
            continue                       # handled as a family below
        lhs = x0.get(a, 0) + tr.get(a, 0)
        # a >= 1 rows carry the +-a twist doubling; a = 0 is its own negative
        target = 2 * _clgr_sum(d) if a >= 1 else _clgr_sum(d)
        tags = ('ram-l' if quad_rec(f, l) == 0 else '') + (',ss' if a == 0 else '')
        rows.append({'a': a, 'd': d, 'lhs': lhs, 'target': target,
                     'delta': lhs - target, 'tags': tags.strip(',')})
    # j = 0, 1728 families: the unit orders O_{-3}, O_{-4} straddle their
    # classes (twists carry distinct traces), so the pairing holds for the
    # family in aggregate, with the straddled order counted once:
    #   sum lhs = 2 [ sum_a clgr_sum(a^2-4p) - (k-1) ],  k = family size
    for f, members in fams.items():
        if not members:
            continue
        # X0 family total: scoped counts plus the j = 0/1728 vertex itself,
        # counted once with its global edge count (as in ecfp_disc_l_isocts)
        qr = quad_rec(f, l)
        if f * l * l in all_ds:
            nu = 1 + qr + (l - qr) // (3 if f == -3 else 2)
        else:
            nu = 1 + qr
        x0_fam = sum(x0.get(a, 0) for a in members) + nu * _clgr_gen(f)
        tr_fam = sum(tr.get(a, 0) for a in members)
        k = len(members)
        target = 2 * (sum(_clgr_sum(a * a - 4 * p) for a in members) - (k - 1))
        rows.append({'a': tuple(members), 'd': f, 'lhs': x0_fam + tr_fam,
                     'target': target, 'delta': x0_fam + tr_fam - target,
                     'tags': f'family j={"0" if f == -3 else "1728"} (k={k})'})
    total = sum(r['lhs'] for r in rows)
    ok_global = (total == 2 * p) and \
                (total == _x0card(p, l) - 2 + _sstr(l, p))
    if verbose:
        print(f'(p={p}, l={l}): per-class + family contributions sum to '
              f'{total} (want {2*p}); consistent with formulas: {ok_global}')
        for r in rows:
            if r['delta'] or r['tags']:
                print(f"  a={str(r['a']):>9s} d={r['d']:6d}  lhs={r['lhs']:4d} "
                      f"target={r['target']:4d} delta={r['delta']:3d}  [{r['tags']}]")
    return {'rows': rows, 'total': total, 'ok_global': ok_global}

def pairing_catalog(N=60, verbose=True):
    """Run pairing_table over all pairs p != l < N (p, l >= 5); summarize
    which (tags -> delta) patterns occur."""
    ps = [q for q in primesBetween(5, N)]
    patterns = {}
    bad_global = []
    for p in ps:
        for l in ps:
            if p == l:
                continue
            t = pairing_table(p, l, verbose=False)
            if not t['ok_global']:
                bad_global.append((p, l, t['total']))
            for r in t['rows']:
                key = (r['tags'], r['delta'])
                patterns[key] = patterns.get(key, 0) + 1
    if verbose:
        print(f'global failures: {len(bad_global)}', bad_global[:10])
        print('observed (tags -> delta) patterns:')
        for (tags, delta), ct in sorted(patterns.items()):
            print(f'  [{tags or "generic"}] delta={delta:3d}   x{ct}')
    return patterns, bad_global

##############################################
# Stage 4: the symmetry grid beyond Atkin l  #
##############################################

_RESULTS_PATH = Path(__file__).parent.parent / 'experiments' / 'pell_duality_results.json'

def _results_load():
    if _RESULTS_PATH.exists():
        with open(_RESULTS_PATH) as f:
            return json.load(f)
    return {}

def _results_save(res):
    _RESULTS_PATH.parent.mkdir(exist_ok=True)
    tmp = _RESULTS_PATH.with_suffix('.tmp')
    with open(tmp, 'w') as f:
        json.dump(res, f, indent=1, sort_keys=True)
    tmp.replace(_RESULTS_PATH)

def stage4_grid(primes, verbose=True):
    """Compute Q_slice/Q_full for all ordered pairs from `primes` (cached)."""
    res = _results_load()
    todo = [(p, l) for p in primes for l in primes
            if p != l and f'{p},{l}' not in res]
    for k, (p, l) in enumerate(todo):
        t0 = time.time()
        out = stage2_Q(p, l)
        out['time'] = round(time.time() - t0, 1)
        res[f'{p},{l}'] = out
        _results_save(res)
        if verbose:
            print(f"  [{k+1}/{len(todo)}] (p={p}, l={l}): "
                  f"Q_slice={out['Q_slice']} Q_full={out['Q_full']} "
                  f"({out['time']}s)")
    return res

def stage4_report(primes):
    """Symmetry defects Q(p,l) - Q(l,p) for both flavors, over `primes`."""
    res = _results_load()
    rows = []
    for i, p in enumerate(primes):
        for l in primes[i + 1 :]:
            a, b = res.get(f'{p},{l}'), res.get(f'{l},{p}')
            if a is None or b is None:
                continue
            rows.append({'p': p, 'l': l,
                         'Qs_pl': a['Q_slice'], 'Qs_lp': b['Q_slice'],
                         'ds': a['Q_slice'] - b['Q_slice'],
                         'Qf_pl': a['Q_full'], 'Qf_lp': b['Q_full'],
                         'df': a['Q_full'] - b['Q_full']})
    return rows


if __name__ == '__main__':
    import sys
    if 'stage1' in sys.argv:
        out = stage1()
        ok = not out['mismatches'] and not out['asymmetries']
        print('STAGE 1:', 'OK — notebook matrix reproduced exactly' if ok else 'DISCREPANCIES FOUND')
    if 'stage2' in sys.argv:
        bad = stage2_calibrate()
        print('STAGE 2:', 'OK — matches stage 1' if not bad else f'{len(bad)} MISMATCHES: {bad}')
