from nt import gcd_list

############
# Matrices #
############

def shape(m: list[list[int]]):
    d2s = list({len(r) for r in m})
    if len(d2s) > 1:
        raise ValueError('Not a matrix')
    return (len(m), d2s[0])

def dot(v1: list, v2: list):
    if len(v1) != len(v2):
        raise ValueError('Vectors must be of equal dimension')
    return sum(x * y for x, y in zip(v1, v2))

def mat_mul(m1, m2):
    d0, d1 = shape(m1)
    d2, d3 = shape(m2)
    if d1 != d2:
        raise ValueError('Incompatible matrices')
    mprod = [[0] * d3 for _ in range(d0)]
    for i, r in enumerate(m1):
        for j in range(d3):
            cj = [r2[j] for r2 in m2]
            mprod[i][j] = dot(r, cj)
    return mprod

def check_square_matrix(m) -> bool:
    s = {len(r) for r in m}
    if len(s) > 1:
        return False
    return min(s) == len(m)

def dotvec(v1, v2):
    if len(v1) != len(v2):
        raise ValueError('Lengths must be equal')
    return sum(x * y for x, y in zip(v1, v2))


class IntegerSquareMatrix:
    def __init__(self, data: list[list[int]]):
        if type(data) != list or type(data[0]) != list:
            raise ValueError('Input must be list of lists')
        if len(data) != len(data[0]):
            raise ValueError('Input must be square matrix')
        self.mat = data
        self.dim = len(data)

    def __str__(self):
        m = self.mat
        s = '['
        for r in m[:-1]:
            s += str(r) + ',\n '
        s += str(m[-1]) + ']'
        return s

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if not isinstance(other, IntegerSquareMatrix):
            return False
        return self.mat == other.mat

    def mvec(self, v: list):
        if len(v) != self.dim:
            raise ValueError('Wrong dimension')
        return [dot(r, v) for r in self.mat]

    def gcd(self):
        return gcd_list([gcd_list(r) for r in self.mat])

    def gcdfac(self):
        d = self.gcd()
        if d < 2:
            return self, d
        return IntegerSquareMatrix([[a // d for a in r] for r in self.mat]), d

    def transpose(self):
        mat = self.mat
        return IntegerSquareMatrix([[r[i] for r in mat] for i in range(len(mat))])

    def trace(self):
        return sum(self.mat[i][i] for i in range(self.dim))

    def det(self):
        if self.dim != 2:
            raise ValueError('Not yet implemented')
        a, b, c, d = self.mat[0][0], self.mat[0][1], self.mat[1][0], self.mat[1][1]
        return a * d - b * c

    def __rmul__(self, n: int):
        if not isinstance(n, int):
            return self * n
        return IntegerSquareMatrix([[x * n for x in r] for r in self.mat])

    def __neg__(self):
        return (-1) * self

    def __add__(self, other: 'IntegerSquareMatrix'):
        if isinstance(other, int):
            other = other * IntegerSquareMatrix([[int(i == j) for j in range(self.dim)]
                                                 for i in range(self.dim)])
        if not isinstance(other, IntegerSquareMatrix):
            raise ValueError('Can only add matrices')
        if other.dim != self.dim:
            raise ValueError('Dimensions must be equal')
        return IntegerSquareMatrix([[self.mat[i][j] + other.mat[i][j]
                                     for j in range(self.dim)]
                                    for i in range(self.dim)])

    def __mul__(self, other: 'IntegerSquareMatrix'):
        if isinstance(other, int):
            other = other * IntegerSquareMatrix([[int(i == j) for j in range(self.dim)]
                                                 for i in range(self.dim)])
        if not isinstance(other, IntegerSquareMatrix):
            raise ValueError('Can only multiply square matrices')
        if other.dim != self.dim:
            raise ValueError('Dimensions must be equal')
        mat2t = other.transpose().mat
        return IntegerSquareMatrix([[dotvec(r1, r2) for r2 in mat2t]
                                    for r1 in self.mat])

    def __pow__(self, n: int):
        if not isinstance(n, int) or n < 0:
            raise ValueError('Can only raise to nonnegative integer powers')
        m = self.dim
        mn = IntegerSquareMatrix([[int(i == j) for j in range(m)] for i in range(m)])
        m2k = self
        while n > 0:
            if n % 2 == 1:
                mn = mn * m2k
            m2k = m2k * m2k
            n //= 2
        return mn

    def __sub__(self, other: 'IntegerSquareMatrix') -> 'IntegerSquareMatrix':
        return self + (-1) * other

    def polyeval(self, coefs_x0_to_xn: list[int]):
        one = self ** 0
        polymat = 0 * one
        for c in coefs_x0_to_xn[::-1]:
            polymat = self * polymat + c * one
        return polymat


def poly2mat(coefs_x0to_xn1: list[int]):
    dim = len(coefs_x0to_xn1)
    mat = [[int(i == j + 1) for j in range(dim)] for i in range(dim)]
    for i in range(dim):
        mat[i][-1] -= coefs_x0to_xn1[i]
    return IntegerSquareMatrix(mat)

def monicpolyremainder(moniccoefs_x0_to_xn, qcoefs):
    mat = poly2mat(moniccoefs_x0_to_xn[:-1])
    polyatmat = mat.polyeval(qcoefs)
    return [r[0] for r in polyatmat.mat]
