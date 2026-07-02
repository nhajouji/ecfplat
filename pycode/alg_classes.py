from typing import Callable
from functools import lru_cache
from nt import gcd_list

                ###################
                # General Classes #
                ###################



class AbGrp:
    def __init__(self,membership:Callable[[tuple],bool],zero:tuple,negation:Callable[[tuple],tuple],addition:Callable[[tuple,tuple],tuple]):
        self.zero_element = zero
        self.negation_map = negation
        self.addition_map = addition
        self.membership_function = membership

    def __contains__(self,elt:tuple):
        if elt == self.zero_element:
            return True
        mem = self.membership_function
        return mem(tuple(elt))
    

    def add_elements(self,t1:tuple,t2:tuple):
        if t1 not in self or t2 not in self:
            raise ValueError('Can only add elements from the same abelian group')
        t0 = self.zero_element
        if t1 == t0:
            return t2
        elif t2 == t0:
            return t1
        neg = self.negation_map
        if neg(t1)==t2:
            return t0
        add = self.addition_map
        return add(t1,t2)

    def negate_element(self, t1:tuple):
        if t1 not in self:
            raise ValueError('Element not in group')
        neg = self.negation_map
        return neg(t1)

    def subtract_elements(self,t1:tuple,t2:tuple):
        t2n = self.negate_element(t2)
        return self.add_elements(t1,t2n)

    def scale_element(self,t1:tuple,n:int):
        if n < 0:
            t1 = self.negate_element(t1)
            n = -n
        t1_n0 = self.zero_element
        t1_n2k = t1
        while n > 0:
            if n % 2 == 1:
                t1_n0 = self.add_elements(t1_n0,t1_n2k)
            n = n//2
            t1_n2k = self.add_elements(t1_n2k,t1_n2k)
        return t1_n0



class AbGrElt:
    def __init__(self,xs:tuple,group:AbGrp):
        if xs not in group:
            raise ValueError('Element not in group')
        self.vec = xs
        self.grp = group
    def __repr__(self):
        return str(self.vec)
    def __add__(self,other:'AbGrElt'):
        group = self.grp
        if not isinstance(other, AbGrElt) or other.vec not in group:
            raise ValueError('Can only add elements of same group')
        t_sum = group.add_elements(self.vec,other.vec)
        return type(self)(t_sum,group)
    def __neg__(self):
        group = self.grp
        t_neg = group.negate_element(self.vec)
        return type(self)(t_neg,group)
    def __sub__(self,other:'AbGrElt'):
        group = self.grp
        if not isinstance(other, AbGrElt) or other.vec not in group:
            raise ValueError('Can only subtract elements of same group')
        t_diff = group.subtract_elements(self.vec,other.vec)
        return type(self)(t_diff,group)
    def __rmul__(self,n:int):
        if type(n)!= int:
            raise ValueError('Can only scale by integers')
        group = self.grp
        tsc = group.scale_element(self.vec,n)
        return type(self)(tsc,group)

class Ring(AbGrp):
    def __init__(
        self,
        membership: Callable[[tuple], bool],
        zero: tuple,
        negation: Callable[[tuple], tuple],
        addition: Callable[[tuple, tuple], tuple],
        one: tuple,
        multiplication: Callable[[tuple, tuple], tuple],
    ):
        super().__init__(membership, zero, negation, addition)
        self.one_element = one
        self.multiplication_map = multiplication

    def multiply_elements(self, t1: tuple, t2: tuple) -> tuple:
        if t1 not in self or t2 not in self:
            raise ValueError('Can only multiply elements from the same ring')
        return self.multiplication_map(t1, t2)


class RingElement(AbGrElt):
    def __init__(self, xs: tuple, ring: Ring):
        super().__init__(xs, ring)

    def __mul__(self, other):
        if isinstance(other, int):
            return self.__rmul__(other)  # scalar scaling from AbGrElt
        if not isinstance(other, RingElement) or other.grp is not self.grp:
            raise ValueError('Can only multiply elements of the same ring')
        product = self.grp.multiply_elements(self.vec, other.vec)
        return type(self)(product, self.grp)

    def __pow__(self, n: int):
        if n < 0:
            raise ValueError('Negative powers not supported (no division)')
        result = type(self)(self.grp.one_element, self.grp)
        base = self
        while n > 0:
            if n % 2 == 1:
                result = result * base
            base = base * base
            n //= 2
        return result


class Field(Ring):
    def __init__(
        self,
        membership: Callable[[tuple], bool],
        zero: tuple,
        negation: Callable[[tuple], tuple],
        addition: Callable[[tuple, tuple], tuple],
        one: tuple,
        multiplication: Callable[[tuple, tuple], tuple],
        inverse: Callable[[tuple], tuple],
    ):
        super().__init__(membership, zero, negation, addition, one, multiplication)
        self.inverse_map = inverse

    def invert_element(self, t: tuple) -> tuple:
        if t == self.zero_element:
            raise ZeroDivisionError('Zero has no inverse')
        if t not in self:
            raise ValueError('Element not in field')
        return self.inverse_map(t)

    def divide_elements(self, t1: tuple, t2: tuple) -> tuple:
        return self.multiply_elements(t1, self.invert_element(t2))


class FieldElement(RingElement):
    def inverse(self) -> 'FieldElement':
        return type(self)(self.grp.invert_element(self.vec), self.grp)

    def __truediv__(self, other: 'FieldElement') -> 'FieldElement':
        if not isinstance(other, FieldElement) or other.grp is not self.grp:
            raise ValueError('Can only divide elements of the same field')
        return self * other.inverse()

    def __pow__(self, n: int) -> 'FieldElement':
        if n < 0:
            return self.inverse() ** (-n)
        return super().__pow__(n)

    def __eq__(self, other) -> bool:
        if not isinstance(other, FieldElement) or other.grp is not self.grp:
            return NotImplemented
        return self.vec == other.vec

    def __hash__(self):
        return hash((id(self.grp), self.vec))

                #######################
                ### Special Classes ###
                #######################

                # Abelian group classes #

def ZnProduct(ns: tuple[int, ...]) -> AbGrp:
    k = len(ns)

    def membership(v: tuple) -> bool:
        return isinstance(v, tuple) and len(v) == k and all(isinstance(x, int) for x in v)

    def negation(v: tuple) -> tuple:
        return tuple(-x % n for x, n in zip(v, ns))

    def addition(v1: tuple, v2: tuple) -> tuple:
        return tuple((x1 + x2) % n for x1, x2, n in zip(v1, v2, ns))

    return AbGrp(membership, tuple(0 for _ in ns), negation, addition)



                        # Rings #

# --- Matrix ring over Z ---

def Mat_n_Z(n: int) -> Ring:
    zero_mat = tuple(tuple(0 for _ in range(n)) for _ in range(n))
    one_mat  = tuple(tuple(1 if i == j else 0 for j in range(n)) for i in range(n))

    def membership(v: tuple) -> bool:
        return (
            isinstance(v, tuple) and len(v) == n
            and all(isinstance(row, tuple) and len(row) == n
                    and all(isinstance(x, int) for x in row)
                    for row in v)
        )

    def negation(v: tuple) -> tuple:
        return tuple(tuple(-x for x in row) for row in v)

    def addition(v1: tuple, v2: tuple) -> tuple:
        return tuple(tuple(x + y for x, y in zip(r1, r2)) for r1, r2 in zip(v1, v2))

    def multiplication(v1: tuple, v2: tuple) -> tuple:
        return tuple(
            tuple(sum(v1[i][k] * v2[k][j] for k in range(n)) for j in range(n))
            for i in range(n)
        )

    return Ring(membership, zero_mat, negation, addition, one_mat, multiplication)


class MatrixElement(RingElement):
    def __init__(self, xs: tuple, ring: Ring):
        super().__init__(xs, ring)
        self.n = len(xs)

    def __mul__(self, other):
        if isinstance(other, int):
            return self.__rmul__(other)
        if not isinstance(other, MatrixElement) or other.grp is not self.grp:
            raise ValueError('Can only multiply elements of the same ring')
        product = self.grp.multiply_elements(self.vec, other.vec)
        return MatrixElement(product, self.grp)

    def __eq__(self, other) -> bool:
        if not isinstance(other, MatrixElement):
            return NotImplemented
        return self.vec == other.vec

    def __repr__(self):
        rows = [' '.join(f'{x:4}' for x in row) for row in self.vec]
        return '\n'.join(rows)

    @property
    def trace(self) -> int:
        return sum(self.vec[i][i] for i in range(self.n))

    @property
    def transpose(self) -> 'MatrixElement':
        t = tuple(tuple(self.vec[j][i] for j in range(self.n)) for i in range(self.n))
        return MatrixElement(t, self.grp)

    @property
    def diagonal(self) -> tuple:
        return tuple(self.vec[i][i] for i in range(self.n))

    @property
    def det(self) -> int:
        if self.n != 2:
            raise NotImplementedError('det only implemented for 2x2 matrices')
        (a, b), (c, d) = self.vec
        return a * d - b * c

    @property
    def adjugate(self) -> 'MatrixElement':
        if self.n != 2:
            raise NotImplementedError('adjugate only implemented for 2x2 matrices')
        (a, b), (c, d) = self.vec
        return MatrixElement(((d, -b), (-c, a)), self.grp)


                        # Fields #

# --- Prime field F_p ---
# Elements are 1-tuples (x,) with 0 <= x < p, so F_p slots into the same
# coefficient-tuple convention as F_{p^n} below.

@lru_cache(maxsize=None)
def GF_p(p: int) -> Field:
    def membership(v: tuple) -> bool:
        return (isinstance(v, tuple) and len(v) == 1
                and isinstance(v[0], int) and 0 <= v[0] < p)

    def negation(v: tuple) -> tuple:
        return ((-v[0]) % p,)

    def addition(v1: tuple, v2: tuple) -> tuple:
        return ((v1[0] + v2[0]) % p,)

    def multiplication(v1: tuple, v2: tuple) -> tuple:
        return ((v1[0] * v2[0]) % p,)

    def inverse(v: tuple) -> tuple:
        return (pow(v[0], p - 2, p),)            # Fermat: x^(p-2) = x^{-1}

    F = Field(membership, (0,), negation, addition, (1,), multiplication, inverse)
    F.char, F.degree, F.order = p, 1, p
    F.int_coercion = lambda c: ((c % p),)
    return F


# --- Extension field F_{p^n} = F_p[x]/(modpoly) ---
# Elements are coefficient n-tuples (a_0, ..., a_{n-1}) for a_0 + a_1 x + ...
# modpoly is a monic irreducible of degree n (low-to-high coefficients); supply a
# Conway polynomial for a canonical, compatible model of F_{p^n}.

def _poly_rem_fp(coefs: list, modpoly: list, p: int) -> list:
    """Exact remainder of a polynomial (low-to-high coefs) modulo a monic modpoly,
    over F_p.  Unlike Polynomial.__mod__ this keeps the true remainder, not its
    monic associate, which is what field reduction needs."""
    r = [c % p for c in coefs]
    m = len(modpoly) - 1                          # deg(modpoly), monic
    while len(r) - 1 >= m:
        while len(r) > 1 and r[-1] == 0:
            r.pop()
        if len(r) - 1 < m:
            break
        lc, shift = r[-1], len(r) - 1 - m
        for i, c in enumerate(modpoly):
            r[i + shift] = (r[i + shift] - lc * c) % p
        r.pop()                                   # leading coefficient is now 0
    return r


def GF_pn(p: int, modpoly: list[int]) -> Field:
    n = len(modpoly) - 1
    if modpoly[-1] != 1:
        raise ValueError('modpoly must be monic')
    order = p ** n
    zero = tuple(0 for _ in range(n))
    one = (1,) + tuple(0 for _ in range(n - 1))

    def membership(v: tuple) -> bool:
        return (isinstance(v, tuple) and len(v) == n
                and all(isinstance(x, int) and 0 <= x < p for x in v))

    def negation(v: tuple) -> tuple:
        return tuple((-x) % p for x in v)

    def addition(v1: tuple, v2: tuple) -> tuple:
        return tuple((x + y) % p for x, y in zip(v1, v2))

    def multiplication(v1: tuple, v2: tuple) -> tuple:
        prod = [0] * (2 * n - 1)
        for i, a in enumerate(v1):
            for j, b in enumerate(v2):
                prod[i + j] += a * b
        rem = _poly_rem_fp(prod, modpoly, p)
        return tuple((rem + [0] * n)[:n])

    def inverse(v: tuple) -> tuple:               # v^(p^n - 2) = v^{-1}
        e, result, base = order - 2, one, v
        while e > 0:
            if e & 1:
                result = multiplication(result, base)
            base = multiplication(base, base)
            e >>= 1
        return result

    K = Field(membership, zero, negation, addition, one, multiplication, inverse)
    K.char, K.degree, K.order = p, n, order
    K.int_coercion = lambda c: ((c % p),) + zero[1:]
    return K


def GF_p2(p: int, modpoly: list[int] = None) -> Field:
    """F_{p^2}.  With no modpoly, uses x^2 - s for the least non-residue s (p odd)."""
    if modpoly is None:
        if p == 2:
            modpoly = [1, 1, 1]                   # x^2 + x + 1
        else:
            s = next(t for t in range(2, p) if pow(t, (p - 1) // 2, p) == p - 1)
            modpoly = [(-s) % p, 0, 1]            # x^2 - s
    return GF_pn(p, modpoly)


# --- Finding an irreducible polynomial of given degree over F_p ---
# Any irreducible degree-n polynomial gives a valid model of F_{p^n}; we do NOT
# need the canonical (Conway) one, because an isogeny's codomain lands back in
# the base field regardless of the model used.  So instead of sourcing Conway
# polynomials we just find one on demand (random search + Rabin's test).

def _poly_mulmod_fp(u: list, v: list, modpoly: list, p: int) -> list:
    prod = [0] * (len(u) + len(v) - 1)
    for i, c1 in enumerate(u):
        for j, c2 in enumerate(v):
            prod[i + j] += c1 * c2
    return _poly_rem_fp(prod, modpoly, p)


def _x_powmod_fp(e: int, modpoly: list, p: int) -> list:
    """x^e mod modpoly over F_p (low-to-high coefficients)."""
    result, base = [1], _poly_rem_fp([0, 1], modpoly, p)
    while e > 0:
        if e & 1:
            result = _poly_mulmod_fp(result, base, modpoly, p)
        base = _poly_mulmod_fp(base, base, modpoly, p)
        e >>= 1
    return result


def irreducible_poly(p: int, n: int, seed: int = 0) -> list[int]:
    """A monic irreducible polynomial of degree n over F_p (low-to-high coefs),
    found by random search + Rabin's irreducibility test."""
    if n == 1:
        return [0, 1]
    import random
    from nt import primefact
    qs = list(primefact(n))                       # distinct primes dividing n
    rng = random.Random(p * 1000003 + n * 101 + seed)
    x = [0, 1]
    while True:
        f = [rng.randrange(p) for _ in range(n)] + [1]            # monic, degree n
        # Rabin: irreducible iff x^(p^n) == x (mod f) and, for each prime q | n,
        # gcd(x^(p^(n/q)) - x, f) == 1.
        if _x_powmod_fp(p ** n, f, p) != _poly_rem_fp(x, f, p):
            continue
        good = True
        for q in qs:
            diff = PolyFp(_x_powmod_fp(p ** (n // q), f, p), p) - PolyFp(x, p)
            if diff.deg < 0 or poly_gcd(diff, PolyFp(f, p)).deg > 0:
                good = False                      # shares a lower-degree factor
                break
        if good:
            return f


@lru_cache(maxsize=None)
def GF_pn_auto(p: int, n: int, seed: int = 0) -> Field:
    """F_{p^n} built on an irreducible polynomial found on demand.

    Memoized per (p, n, seed): irreducible_poly is a deterministic (seeded) search,
    so reusing the field is sound and skips the Rabin search + field rebuild that
    otherwise ran on every Velu isogeny call.  Reusing one field object per (p, n)
    also keeps FieldElement equality/hashing consistent across calls."""
    return GF_pn(p, irreducible_poly(p, n, seed))


                # Coefficient rings for polynomial work #

# Z and C as Ring/Field instances, so PolyRing below can be built over them with
# the same element-tuple convention as GF_p / GF_pn.  CC is floating point: it is
# meant for numeric estimates (e.g. coefficient bounds from approximate CM
# j-invariants), not exact arithmetic.

def _make_ZZ() -> Ring:
    def membership(v):
        return isinstance(v, tuple) and len(v) == 1 and isinstance(v[0], int)
    R = Ring(membership, (0,), lambda v: (-v[0],), lambda v1, v2: (v1[0] + v2[0],),
             (1,), lambda v1, v2: (v1[0] * v2[0],))
    R.int_coercion = lambda c: (c,)
    return R

ZZ = _make_ZZ()


def _make_CC() -> Field:
    def membership(v):
        return isinstance(v, tuple) and len(v) == 1 and isinstance(v[0], (int, float, complex))
    F = Field(membership, (0,), lambda v: (-v[0],), lambda v1, v2: (v1[0] + v2[0],),
              (1,), lambda v1, v2: (v1[0] * v2[0],), lambda v: (1 / v[0],))
    F.int_coercion = lambda c: (c,)
    return F

CC = _make_CC()


                    ####################
                    # Polynomial rings #
                    ####################

# R[x] over any Ring R from this module (ZZ, GF_p, GF_pn, CC, ...).  Elements are
# tuples of base-ring element tuples, low-to-high degree, with trailing zeros
# stripped; the zero polynomial is the empty tuple ().  Build rings through
# poly_ring(base) so equal bases share one PolyRing object.

class PolyRing(Ring):
    def __init__(self, base: Ring):
        self.base = base
        bz, badd, bneg, bmul = (base.zero_element, base.addition_map,
                                base.negation_map, base.multiplication_map)

        def strip(coefs: tuple) -> tuple:
            n = len(coefs)
            while n > 0 and coefs[n - 1] == bz:
                n -= 1
            return tuple(coefs[:n])

        def membership(v):
            return isinstance(v, tuple) and all(c in base for c in v)

        def negation(v):
            return tuple(bneg(c) for c in v)

        def addition(v1, v2):
            if len(v1) < len(v2):
                v1, v2 = v2, v1
            out = list(v1)
            for i, c in enumerate(v2):
                out[i] = badd(out[i], c)
            return strip(tuple(out))

        def multiplication(v1, v2):
            if not v1 or not v2:
                return ()
            out = [bz] * (len(v1) + len(v2) - 1)
            for i, c1 in enumerate(v1):
                if c1 == bz:
                    continue
                for j, c2 in enumerate(v2):
                    out[i + j] = badd(out[i + j], bmul(c1, c2))
            return strip(tuple(out))

        self._strip = strip
        super().__init__(membership, (), negation, addition, (base.one_element,), multiplication)

    # --- element construction ---

    def coerce_coef(self, c) -> tuple:
        """Turn an int, a base-ring tuple, or a base-ring element into a base tuple."""
        if isinstance(c, AbGrElt):
            c = c.vec
        if isinstance(c, tuple):
            if c in self.base:
                return c
            raise ValueError(f'{c} is not in the coefficient ring')
        if isinstance(c, (int, float, complex)) and not isinstance(c, bool):
            coerce = getattr(self.base, 'int_coercion', None)
            if coerce is not None:
                return coerce(c)
            return self.base.scale_element(self.base.one_element, c)
        raise ValueError(f'Cannot coerce {c} into the coefficient ring')

    def poly(self, coefs) -> 'Poly':
        """Polynomial from low-to-high coefficients (ints or base-ring elements)."""
        return Poly(self._strip(tuple(self.coerce_coef(c) for c in coefs)), self)

    @property
    def zero(self) -> 'Poly':
        return Poly((), self)

    @property
    def one(self) -> 'Poly':
        return Poly(self.one_element, self)

    @property
    def x(self) -> 'Poly':
        return Poly((self.base.zero_element, self.base.one_element), self)

    def from_roots(self, roots) -> 'Poly':
        """Monic product (x - r_1)...(x - r_n), multiplied as a balanced tree."""
        bneg, bone = self.base.negation_map, self.base.one_element
        layer = [(bneg(self.coerce_coef(r)), bone) for r in roots]
        if not layer:
            return self.one
        while len(layer) > 1:
            nxt = [self.multiply_elements(layer[i], layer[i + 1])
                   for i in range(0, len(layer) - 1, 2)]
            if len(layer) % 2:
                nxt.append(layer[-1])
            layer = nxt
        return Poly(layer[0], self)

    # --- division-based operations (coefficients must form a field) ---

    def divmod_polys(self, f: 'Poly', g: 'Poly') -> tuple['Poly', 'Poly']:
        if not isinstance(self.base, Field):
            raise ValueError('Polynomial division needs field coefficients')
        if not g.vec:
            raise ZeroDivisionError('Division by the zero polynomial')
        badd, bneg, bmul = (self.base.addition_map, self.base.negation_map,
                            self.base.multiplication_map)
        dg = len(g.vec) - 1
        lc_inv = self.base.invert_element(g.vec[-1])
        rem = list(f.vec)
        quo = [self.base.zero_element] * max(len(rem) - dg, 0)
        while len(rem) - 1 >= dg:
            if rem[-1] == self.base.zero_element:
                rem.pop()
                continue
            q = bmul(rem[-1], lc_inv)
            shift = len(rem) - 1 - dg
            quo[shift] = q
            for i, c in enumerate(g.vec):
                rem[i + shift] = badd(rem[i + shift], bneg(bmul(q, c)))
            rem.pop()
        return Poly(self._strip(tuple(quo)), self), Poly(self._strip(tuple(rem)), self)

    def gcd_polys(self, f: 'Poly', g: 'Poly') -> 'Poly':
        """Monic gcd (field coefficients)."""
        while g.vec:
            f, g = g, self.divmod_polys(f, g)[1]
        return f.monic() if f.vec else f

    def resultant(self, f: 'Poly', g: 'Poly'):
        """Res(f, g) as a base-field element tuple, via the Euclidean algorithm."""
        if not isinstance(self.base, Field):
            raise ValueError('resultant needs field coefficients')
        base = self.base
        if not f.vec or not g.vec:
            return base.zero_element if max(f.deg, g.deg) > 0 else base.one_element
        res = base.one_element
        while True:
            if not g.vec:
                return base.zero_element if f.deg > 0 else res
            if g.deg == 0:
                return base.multiply_elements(res, _base_pow(base, g.vec[0], f.deg))
            r = self.divmod_polys(f, g)[1]
            dr = r.deg if r.vec else 0
            res = base.multiply_elements(res, _base_pow(base, g.vec[-1], f.deg - dr))
            if (f.deg * g.deg) % 2:
                res = base.negation_map(res)
            f, g = g, r

    def derivative(self, f: 'Poly') -> 'Poly':
        sc = self.base.scale_element
        dcs = tuple(sc(c, i) for i, c in enumerate(f.vec))[1:]
        return Poly(self._strip(dcs), self)


def _base_pow(base: Ring, t: tuple, n: int) -> tuple:
    out, sq = base.one_element, t
    while n > 0:
        if n & 1:
            out = base.multiply_elements(out, sq)
        sq = base.multiply_elements(sq, sq)
        n >>= 1
    return out


_POLY_RINGS: dict = {}

def poly_ring(base: Ring) -> PolyRing:
    """The polynomial ring over base, one shared object per base ring."""
    R = _POLY_RINGS.get(id(base))
    if R is None:
        R = _POLY_RINGS[id(base)] = PolyRing(base)
    return R


class Poly(RingElement):
    @property
    def ring(self) -> PolyRing:
        return self.grp

    @property
    def deg(self) -> int:
        return len(self.vec) - 1          # zero polynomial has degree -1

    @property
    def lc(self) -> tuple:
        return self.vec[-1] if self.vec else self.grp.base.zero_element

    def coef(self, i) -> tuple:
        """Coefficient of x^i as a base-ring tuple (zero beyond the degree)."""
        return self.vec[i] if 0 <= i <= self.deg else self.grp.base.zero_element

    def int_coefs(self) -> list:
        """Low-to-high coefficients as plain scalars (for 1-tuple bases: ZZ, GF_p, CC)."""
        return [c[0] for c in self.vec]

    def __call__(self, x):
        """Evaluate by Horner. x may be an int, a base element, or a Poly (composition)."""
        R = self.grp
        if isinstance(x, Poly):
            out = R.zero
            for c in reversed(self.vec):
                out = out * x + Poly((c,), R)
            return out
        x = R.coerce_coef(x)
        base = R.base
        out = base.zero_element
        for c in reversed(self.vec):
            out = base.addition_map(base.multiplication_map(out, x), c)
        return out

    def monic(self) -> 'Poly':
        if not self.vec or self.lc == self.grp.base.one_element:
            return self
        lc_inv = self.grp.base.invert_element(self.lc)
        bmul = self.grp.base.multiplication_map
        return Poly(tuple(bmul(c, lc_inv) for c in self.vec), self.grp)

    def derivative(self) -> 'Poly':
        return self.grp.derivative(self)

    def __divmod__(self, other: 'Poly'):
        return self.grp.divmod_polys(self, other)

    def __floordiv__(self, other: 'Poly') -> 'Poly':
        return self.grp.divmod_polys(self, other)[0]

    def __mod__(self, other: 'Poly') -> 'Poly':
        return self.grp.divmod_polys(self, other)[1]

    def __eq__(self, other) -> bool:
        if not isinstance(other, Poly) or other.grp.base is not self.grp.base:
            return NotImplemented
        return self.vec == other.vec

    def __hash__(self):
        return hash((id(self.grp.base), self.vec))

    def __repr__(self):
        if not self.vec:
            return '0'
        terms = []
        for i in range(self.deg, -1, -1):
            c = self.vec[i]
            if c == self.grp.base.zero_element:
                continue
            cs = str(c[0]) if len(c) == 1 else str(c)
            xi = '' if i == 0 else ('x' if i == 1 else f'x^{i}')
            if xi and cs == '1':
                terms.append(xi)
            elif xi and cs == '-1':
                terms.append(f'-{xi}')
            else:
                terms.append(f'{cs} {xi}'.strip())
        s = terms[0]
        for t in terms[1:]:
            s += f' - {t[1:]}' if t.startswith('-') else f' + {t}'
        return s


            # Reduction, lifting, and CRT for integer polynomials #

def poly_mod_p(f: Poly, p: int) -> Poly:
    """Reduce a ZZ-polynomial mod p, landing in GF_p(p)[x]."""
    return poly_ring(GF_p(p)).poly(f.int_coefs())


def poly_lift(f: Poly, balanced: bool = True) -> Poly:
    """Lift a GF_p-polynomial to ZZ[x]; balanced puts coefficients in (-p/2, p/2]."""
    p = f.grp.base.char
    cs = f.int_coefs()
    if balanced:
        cs = [c - p if 2 * c > p else c for c in cs]
    return poly_ring(ZZ).poly(cs)


def poly_crt(residues: list[tuple]) -> tuple[Poly, int]:
    """CRT a list of (polynomial mod p, p) pairs into (F, M) with F in ZZ[x], M = prod p.

    Each polynomial may be a Poly over GF_p or a low-to-high list of int coefficients.
    Coefficients of F are the balanced representatives in (-M/2, M/2], which are the
    true integer coefficients whenever M exceeds twice their largest absolute value."""
    from nt import crt_list
    pairs = []
    for f, p in residues:
        cs = f.int_coefs() if isinstance(f, Poly) else list(f)
        pairs.append((cs, p))
    n = max(len(cs) for cs, _ in pairs)
    coefs = []
    for i in range(n):
        ci, M = crt_list([(cs[i] if i < len(cs) else 0, p) for cs, p in pairs])
        coefs.append(ci)
    return poly_ring(ZZ).poly(coefs), M


                    ######################
                    # Legacy polynomials #
                    ######################
# Superseded by PolyRing/Poly above; still used by modularpolynomials.py and
# irreducible_poly.  New code should use poly_ring(...).

def remove_lead_zeros(l, char=0):
    if char != 0:
        l = [c % char for c in l]
    if len(l) == 0:
        return [0]
    while len(l) > 1 and l[-1] == 0:
        l = l[:-1]
    return l

# TODO: Rewrite this so that it follows previous pattern
class Polynomial:
    def __init__(self, coefs_0_to_n: list[int], char=0):
        self.coefs = remove_lead_zeros(coefs_0_to_n, char)
        self.char = char
        self.lc = (self.coefs)[-1]
        self.deg = len(self.coefs) - 1 - int(self.lc == 0)
        self.const = (self.coefs)[0]

    def __str__(self):
        coefs = self.coefs
        if len(coefs) == 1:
            return str(coefs[0])
        s = ''
        for i, c in enumerate(coefs):
            if c != 0:
                if len(s) > 0:
                    if s[0] == '-':
                        s = f'{c} x^{i}' + s
                    else:
                        s = f'{c} x^{i}+' + s
                else:
                    s = f'{c} x^{i}'
        return s if s else '0'

    def __repr__(self):
        return str(self)

    def __add__(self, other: 'Polynomial'):
        if isinstance(other, int):
            other = Polynomial([other], self.char)
        if not isinstance(other, Polynomial) or self.char != other.char:
            raise ValueError('Can only add polynomials of equal char')
        cs1 = self.coefs
        cs2 = other.coefs
        l1, l2 = len(cs1), len(cs2)
        cs1 += max(l2 - l1, 0) * [0]
        cs2 += max(l1 - l2, 0) * [0]
        cssum = remove_lead_zeros([cs1[i] + cs2[i] for i in range(max(l1, l2))], self.char)
        return Polynomial(cssum, self.char)

    def __rmul__(self, n: int):
        csn = [n * c for c in self.coefs]
        if self.char != 0:
            csn = [c % self.char for c in csn]
        return Polynomial(csn, self.char)

    def __sub__(self, other: 'Polynomial'):
        return self + (-1) * other

    def __mul__(self, other: 'Polynomial'):
        if not isinstance(other, (Polynomial, int)):
            raise ValueError('Can only multiply by polynomials or integers')
        if isinstance(other, int):
            return Polynomial([c * other for c in self.coefs], self.char)
        if self.char != other.char:
            raise ValueError('Characteristics must be equal')
        cs1, cs2 = self.coefs, other.coefs
        cs12 = [0] * (len(cs1) + len(cs2) - 1)
        for i, c1 in enumerate(cs1):
            for j, c2 in enumerate(cs2):
                cs12[i + j] += c1 * c2
        if self.char > 0:
            cs12 = [c % self.char for c in cs12]
        return Polynomial(cs12, self.char)

    def __pow__(self, n):
        if n < 0:
            raise ValueError('Exponent must be nonnegative')
        xn = Polynomial([1], self.char)
        x2n = self
        while n > 0:
            if n % 2 == 1:
                xn = xn * x2n
            n = n // 2
            x2n = x2n * x2n
        return xn

    def monic_associate(self):
        if self.lc < 0:
            self = (-1) * self
        if self.deg == -1:
            return self
        elif self.lc == 1:
            return self
        elif self.char > 0 and self.lc != 0:
            lcinv = pow(self.lc, -1, self.char)
            return lcinv * self
        else:
            g = gcd_list(self.coefs)
            return Polynomial([c // g for c in self.coefs])

    def __mod__(self, den: 'Polynomial') -> 'Polynomial':
        if not isinstance(den, Polynomial) or den.char != self.char:
            raise ValueError('Can only mod polynomials of equal char')
        if den.deg < 0:
            raise ZeroDivisionError('Denominator must be nonzero')
        if den.deg == 0:
            return 0 * self
        den0 = den.monic_associate()
        rem0 = self.monic_associate()
        px = Polynomial([0, 1], self.char)
        while rem0.deg >= den0.deg:
            lcr, lcd = rem0.lc, den0.lc
            dd = rem0.deg - den0.deg
            rem_new = (lcd * rem0 - (px ** dd) * lcr * den0).monic_associate()
            if rem_new.deg < rem0.deg:
                rem0 = rem_new
            else:
                return rem0
        return rem0

    def mod(self, p: int):
        if self.char > 0 and self.char != p:
            raise ValueError('Already in positive characteristic')
        return PolyFp(self.coefs, char=p)

    def eval(self, x):
        p = self.char
        coefs = self.coefs
        evx = 0 * x
        one = x ** 0
        if p > 0:
            for c in coefs[::-1]:
                evx = (evx * x + c * one) % p
        else:
            for c in coefs[::-1]:
                evx = evx * x + c * one
        return evx

    def dx(self):
        if len(self.coefs) == 1:
            return 0 * self
        dcs = [i * c for i, c in enumerate(self.coefs)][1:]
        return Polynomial(dcs, char=self.char)


def poly_gcd(poly1, poly2):
    r = poly2 % poly1
    while r.deg >= 0:
        poly2 = poly1
        poly1 = r
        r = poly2 % poly1
    return poly1


class PolyFp(Polynomial):
    def __init__(self, coefs, char):
        super().__init__(coefs, char)

    def fp_factor(self):
        p = self.char
        xp = PolyFp([0, 1], p)
        fpp = xp ** p - xp
        return poly_gcd(self, fpp)

    def no_fp_roots(self):
        return self.fp_factor().deg

    def find_roots_BrFo(self):
        p = self.char
        return [x for x in range(p) if self.eval(x) == 0]

