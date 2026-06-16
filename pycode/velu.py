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
