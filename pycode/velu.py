"""
Velu's formulas: l-isogeny codomains from a kernel, over any field K.

Elliptic curve point arithmetic and Velu's formulas are written over an arbitrary
field built with the Field classes in alg_classes, so the same code runs over F_p
and over extensions F_{p^k}.  This is the route to l-isogeny neighbour data for
primes l where no modular polynomial is available (and, eventually, for
supersingular classes).

Curves are short Weierstrass y^2 = x^3 + a*x + b with a, b given as FieldElements.
Points are either None (the point at infinity) or an (x, y) pair of FieldElements.
"""

from alg_classes import FieldElement


def _zero(K):
    return FieldElement(K.zero_element, K)


#########################
# EC point arithmetic   #
#########################

def ec_neg(P):
    if P is None:
        return None
    x, y = P
    return (x, -y)


def ec_add(P, Q, a):
    """Group law on y^2 = x^3 + a*x + b (b is not needed for addition)."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    K = a.grp
    if x1 == x2:
        if (y1 + y2) == _zero(K):           # Q == -P (covers 2-torsion: y1 == 0)
            return None
        lam = (3 * (x1 * x1) + a) / (2 * y1)        # doubling
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - x1 - x2
    y3 = lam * (x1 - x3) - y1
    return (x3, y3)


def ec_double(P, a):
    return ec_add(P, P, a)


def ec_mul(n, P, a):
    """Scalar multiple n*P by double-and-add."""
    if n < 0:
        return ec_mul(-n, ec_neg(P), a)
    R, base = None, P
    while n > 0:
        if n & 1:
            R = ec_add(R, base, a)
        base = ec_double(base, a)
        n >>= 1
    return R


def point_order(P, a, bound):
    """Smallest k > 0 with k*P == O, searching up to bound (None if not found)."""
    R = P
    for k in range(1, bound + 1):
        if R is None:
            return k
        R = ec_add(R, P, a)
    return None


#########################
# Velu's formulas       #
#########################

def kernel_reps(P, l, a):
    """Representatives of (<P> \\ {O}) / {+-1} for a kernel generator P of order l."""
    if l == 2:
        return [P]                           # the single 2-torsion point
    reps, Q = [], P
    for _ in range((l - 1) // 2):
        reps.append(Q)
        Q = ec_add(Q, P, a)
    return reps


def velu_codomain(a, b, reps):
    """Codomain (A, B) of the isogeny with the given kernel representatives.

    Velu: for each Q = (xQ, yQ) in the representative set,
        gxQ = 3 xQ^2 + a
        tQ  = gxQ            and uQ = 0          if Q is 2-torsion (yQ = 0)
        tQ  = 2 gxQ          and uQ = 4 yQ^2     otherwise
    then A = a - 5*sum(tQ),  B = b - 7*sum(uQ + xQ*tQ)."""
    K = a.grp
    zero = _zero(K)
    t_sum, w_sum = zero, zero
    for (xq, yq) in reps:
        gxq = 3 * (xq * xq) + a
        if yq == zero:
            tq, uq = gxq, zero
        else:
            tq, uq = 2 * gxq, 4 * (yq * yq)
        t_sum = t_sum + tq
        w_sum = w_sum + (uq + xq * tq)
    A = a - 5 * t_sum
    B = b - 7 * w_sum
    return A, B


def j_invariant(a, b):
    """j-invariant of y^2 = x^3 + a*x + b as a FieldElement."""
    num = 4 * (a ** 3)
    return 1728 * num / (num + 27 * (b ** 2))


def velu_isogeny(a, b, l, kgen):
    """Codomain (A, B) of the degree-l isogeny with kernel <kgen>.

    kgen is a kernel generator (a point of order l).  Pure field arithmetic over
    kgen's field, so it runs over F_p or any extension."""
    return velu_codomain(a, b, kernel_reps(kgen, l, a))


#####################################
# Finding the l-torsion kernel      #
#####################################
# To run Velu we need a generator of the kernel eigenline.  We work in F_{p^k}
# (k from frob_ext_degrees), take a random point, multiply by the cofactor to
# land in the l-torsion, then project onto the lambda-eigenline with Frobenius:
#     Q = pi(P) - mu*P,    pi(x,y) = (x^p, y^p),   mu = p/lambda mod l.
# Q (when != O) generates the kernel for the +x_l direction; the codomain j of a
# horizontal isogeny lands back in F_p.

import random as _random


def _rand_elt(K, rng):
    return FieldElement(tuple(rng.randrange(K.char) for _ in range(K.degree)), K)


def embed_fp(c, K):
    """Embed an F_p element (int) into K = F_{p^k}."""
    return FieldElement((c % K.char,) + tuple(0 for _ in range(K.degree - 1)), K)


def field_sqrt(c, rng=None):
    """A square root of c in its field, or None if c is a non-square (Tonelli-Shanks)."""
    K = c.grp
    q = K.order
    zero, one = _zero(K), FieldElement(K.one_element, K)
    if c == zero:
        return zero
    if c ** ((q - 1) // 2) != one:
        return None
    s, t = 0, q - 1
    while t % 2 == 0:
        s, t = s + 1, t // 2
    if s == 1:                                   # q == 3 mod 4
        return c ** ((q + 1) // 4)
    rng = rng or _random.Random(0)
    z = _rand_elt(K, rng)                         # a quadratic non-residue
    while z == zero or z ** ((q - 1) // 2) == one:
        z = _rand_elt(K, rng)
    M, cc, R, tt = s, z ** t, c ** ((t + 1) // 2), c ** t
    while tt != one:
        i, temp = 1, tt * tt
        while temp != one:
            temp, i = temp * temp, i + 1
        b = cc ** (2 ** (M - i - 1))
        M, cc, R, tt = i, b * b, R * b, tt * b * b
    return R


def random_point(a, b, rng):
    """A random affine point on y^2 = x^3 + a*x + b over a.grp."""
    while True:
        x = _rand_elt(a.grp, rng)
        y = field_sqrt(x * x * x + a * x + b, rng)
        if y is not None:
            return (x, y)


def frobenius(P, p):
    """The p-power Frobenius pi(x, y) = (x^p, y^p)."""
    if P is None:
        return None
    x, y = P
    return (x ** p, y ** p)


def curve_cardinality(a, p, k):
    """#E(F_{p^k}) = p^k + 1 - a_k, where a_0=2, a_1=a, a_i = a*a_{i-1} - p*a_{i-2}."""
    t0, t1 = 2, a
    for _ in range(2, k + 1):
        t0, t1 = t1, a * t1 - p * t0
    return p ** k + 1 - t1


def kernel_gen_eigenline(a, b, l, lam, mu, p, k, trace, rng, tries=300):
    """A generator of the lambda-eigenline l-torsion kernel of y^2=x^3+ax+b over
    F_{p^k} (a, b FieldElements).  trace is the F_p trace of Frobenius."""
    N = curve_cardinality(trace, p, k)
    m = N                                       # remove the full l-part of N, not just
    while m % l == 0:                           # one factor (the l-Sylow can be (Z/l)^2,
        m //= l                                 # where N/l would kill all l-torsion)
    mu_red = mu % l
    for _ in range(tries):
        T = ec_mul(m, random_point(a, b, rng), a)        # into the l-Sylow
        if T is None:
            continue
        while ec_mul(l, T, a) is not None:               # reduce to an order-l point
            T = ec_mul(l, T, a)
        Q = ec_add(frobenius(T, p), ec_neg(ec_mul(mu_red, T, a)), a)   # pi(T) - mu*T
        if Q is not None and point_order(Q, a, l) == l:
            return Q
    return None


def velu_l_isog_codomain(f, g, l, p, trace, direction=0, max_degree=None, seed=0):
    """Codomain j of a horizontal l-isogeny of y^2 = x^3 + f*x + g over F_p.

    trace is the F_p trace of Frobenius.  direction selects the eigenvalue (0 or 1
    for split).  Returns a dict with 'status' = 'ok' (and 'j', 'k', 'lam'),
    'inert'/'ramified' (no horizontal l-isogeny), or 'extension_too_large'."""
    from nt import frob_ext_degrees
    from alg_classes import GF_p, GF_pn_auto
    info = frob_ext_degrees(trace, p, l)
    if info['kind'] != 'split':
        return {'status': info['kind']}                  # no two horizontal directions
    eigs, degs = info['eigenvalues'], info['degrees']
    lam, k = eigs[direction], degs[direction]
    mu = (p * pow(lam, -1, l)) % l                        # the other eigenvalue
    if max_degree is not None and k > max_degree:
        return {'status': 'extension_too_large', 'k': k}
    K = GF_p(p) if k == 1 else GF_pn_auto(p, k)
    a, b = embed_fp(f, K), embed_fp(g, K)
    rng = _random.Random((p * 7919 + f) * 7919 + g + l * 131 + seed)
    Q = kernel_gen_eigenline(a, b, l, lam, mu, p, k, trace, rng)
    if Q is None:
        return {'status': 'kernel_not_found', 'k': k}
    A, B = velu_isogeny(a, b, l, Q)
    j = j_invariant(A, B)
    jval = j.vec[0] if all(c == 0 for c in j.vec[1:]) else j.vec     # horizontal -> in F_p
    return {'status': 'ok', 'j': jval, 'k': k, 'lam': lam, 'mu': mu}


def velu_nbr_data_ord(ap, l, js=None, max_degree=None):
    """Horizontal l-isogeny neighbour graph {j: [neighbour j's]} for the ordinary
    isogeny class (a, p), via Velu.  Mirrors ecfp_nbr_data_ord_X1 (the Atkin path)
    but works for any split prime l, not just the 15 with a modular polynomial.

    For each class j, takes a model of trace a (trfr_to_models) and computes both
    horizontal directions; only codomains that land back in the class are kept.
    max_degree caps the extension F_{p^k} (directions needing a larger k are skipped)."""
    a, p = ap
    from ecfp import trfr_to_js, trfr_to_models
    if js is None:
        js = trfr_to_js(a, p)
    models = trfr_to_models(a, p)              # {j: (f,g) with trace exactly a}
    jset = set(js)
    out = {}
    for j in js:
        fg = models.get(j) if isinstance(models, dict) else None
        nbrs = []
        if fg is not None:
            f, g = fg
            for d in (0, 1):
                res = velu_l_isog_codomain(f, g, l, p, a, direction=d, max_degree=max_degree)
                if res['status'] == 'ok' and res['j'] in jset:
                    nbrs.append(res['j'])
        out[j] = list(dict.fromkeys(nbrs))
    return out
