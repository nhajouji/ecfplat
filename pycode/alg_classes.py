from typing import Callable

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

