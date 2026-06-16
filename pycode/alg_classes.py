from typing import Callable
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

    return Field(membership, (0,), negation, addition, (1,), multiplication, inverse)


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

    return Field(membership, zero, negation, addition, one, multiplication, inverse)


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


def GF_pn_auto(p: int, n: int, seed: int = 0) -> Field:
    """F_{p^n} built on an irreducible polynomial found on demand."""
    return GF_pn(p, irreducible_poly(p, n, seed))


                    ################
                    # Polynomials  #
                    ################

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

