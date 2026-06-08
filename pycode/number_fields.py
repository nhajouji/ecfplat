from nt import get_rou_mod, gcd
from ringclasses import poly2mat

#####################
# Quadratic Integers #
#####################

class QuadraticInteger:
    def __init__(self, elementcoefs: tuple[int], ringpoly: tuple[int]):
        self.poly = ringpoly
        self.matrix = poly2mat(ringpoly).polyeval(elementcoefs)
        self.vec = (self.matrix.transpose().mat)[0]

    def __str__(self):
        x, y = self.vec
        c0, c1 = self.poly
        return f'({x})+({y} X) mod (X^2 + {c1} X + {c0})'

    def __repr__(self):
        return str(self)

    def __add__(self, other):
        if not isinstance(other, QuadraticInteger):
            if isinstance(other, int):
                return self + QuadraticInteger((other, 0), self.poly)
            raise TypeError('Incompatible types')
        if self.poly != other.poly:
            raise TypeError('Use same polynomial')
        vec1, vec2 = self.vec, other.vec
        return QuadraticInteger((vec1[0] + vec2[0], vec1[1] + vec2[1]), self.poly)

    def __mul__(self, other):
        if not isinstance(other, QuadraticInteger):
            if isinstance(other, int):
                return self * QuadraticInteger((other, 0), self.poly)
            raise TypeError('Incompatible types')
        if self.poly != other.poly:
            raise TypeError('Use same polynomial')
        matprod = (self.matrix) * (other.matrix)
        vecprod = (matprod.transpose().mat)[0]
        return QuadraticInteger(vecprod, self.poly)

    def __rmul__(self, n: int):
        x, y = self.vec
        return QuadraticInteger((n * x, n * y), self.poly)

    def __sub__(self, other):
        return self + (-1) * other

    def trace(self):
        return self.matrix.trace()

    def conjugate(self):
        t = self.trace()
        return QuadraticInteger((t, 0), self.poly) - self

    def norm(self):
        return ((self * self.conjugate()).vec)[0]

    def dot(self, other):
        if not isinstance(other, QuadraticInteger) or self.poly != other.poly:
            raise TypeError('Use same polynomial')
        return (self * other.conjugate()).trace()


#####################
# Gaussian Integers #
#####################

def _smallest_rem(a, p):
    a0 = a % p
    a1 = a0 - p
    return a1 if abs(a1) < a0 else a0

def _smallest_q(a, p):
    return (a - _smallest_rem(a, p)) // p

def _improve_approx_gauss(xy, p):
    x, y = xy
    k = p // abs(y)
    x1, x2 = (x * k) % p, (x * (k + 1)) % p
    y1, y2 = p - k * abs(y), (k + 1) * abs(y) - p
    v1, v2 = (x1, y1), (x2, y2)
    n1, n2 = x1**2 + y1**2, x2**2 + y2**2
    return v1 if (n1 == p or n1 <= n2) else v2


class GaussianInteger:
    def __init__(self, ab: tuple[int, int]):
        self.vec = ab
        self.re = ab[0]
        self.im = ab[1]
        self.qi = QuadraticInteger(ab, (1, 0))
        self.tr = self.qi.trace()
        self.nrm = self.qi.norm()

    def __str__(self):
        if self.nrm == 0:
            return '0'
        a, b = self.vec
        if b == 0:
            return str(a)
        if a == 0:
            return str(b) + 'i'
        return str(a) + ('+ i' if b > 0 else '- i') + str(abs(b))

    def __repr__(self):
        return str(self)

    def __add__(self, other):
        x1, y1 = self.vec
        x2, y2 = other.vec
        return GaussianInteger((x1 + x2, y1 + y2))

    def __sub__(self, other):
        x1, y1 = self.vec
        x2, y2 = other.vec
        return GaussianInteger((x1 - x2, y1 - y2))

    def __rmul__(self, n: int):
        x1, y1 = self.vec
        return GaussianInteger((n * x1, n * y1))

    def __mul__(self, other):
        x1, y1 = self.vec
        x2, y2 = other.vec
        return GaussianInteger((x1 * x2 - y1 * y2, x1 * y2 + x2 * y1))

    def conj(self):
        x, y = self.vec
        return GaussianInteger((x, -y))

    def tracedot(self, other):
        return ((self * other.conj()).tr) // 2

    def __floordiv__(self, other):
        if not isinstance(other, GaussianInteger):
            raise ValueError('Can only divide by other Gaussian integers')
        if other.nrm == 0:
            raise ZeroDivisionError('Denominator must have nonzero norm')
        v1 = other
        v2 = GaussianInteger((0, 1)) * v1
        c1 = _smallest_q(v1.tracedot(self), v1.tracedot(v1))
        c2 = _smallest_q(v2.tracedot(self), v2.tracedot(v2))
        return GaussianInteger((c1, c2))

    def __mod__(self, other):
        return self - (self // other) * other


def gcd_gauss(z1: GaussianInteger, z2: GaussianInteger) -> GaussianInteger:
    if z1.nrm == 0:
        return z2
    while (z2 % z1).nrm != 0:
        r = z2 % z1
        z2 = z1
        z1 = r
    return z1


def sum_of_squares_rep(p):
    if p % 4 == 3:
        return None
    if p == 2:
        return GaussianInteger((1, 1))
    r0 = min(get_rou_mod(4, p), p - get_rou_mod(4, p))
    z1 = GaussianInteger((1, r0))
    if z1.nrm == p:
        return z1
    z0 = GaussianInteger((1, p - r0))
    while z1.nrm != p and z1.nrm < z0.nrm:
        if gcd(z0.nrm, z1.nrm) < z1.nrm:
            z01 = gcd_gauss(z0, z1)
            z01c = gcd_gauss(z0, z1.conj())
            if z01.nrm % p == 0:
                z1 = z01 if z01.nrm == p else z01
                if z01.nrm == p:
                    return z01
                z1 = z01
            elif z01c.nrm % p == 0:
                if z01c.nrm == p:
                    return z01c
                z1 = z01c
        if z1.nrm == p:
            return z1
        z2 = GaussianInteger(_improve_approx_gauss(z1.vec, p))
        z0 = z1
        z1 = z2
    return z1
