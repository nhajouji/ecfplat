"""Binary quadratic forms and the SL2(Z) / Gamma_0(l) action.

A form is a tuple (a, b, c) of integers, standing for a x^2 + b xy + c y^2 --
equivalently the point tau = (-b + sqrt(D))/2a of the upper half-plane and the
CM lattice <1, tau>.  Provides: reduction into the fundamental domain
(qf_mod_gamma / _qf_reduce), enumeration of the class group of a discriminant
(get_qfs_strict / get_qfs_all / qfs_ordered_by_cond), l-isogeny neighbours
(qf_isogs_hor = horizontal, qf_parents = ascending; cached), isogeny cycles,
Gamma_0(l) coset orbits, the Fricke involution, and endomorphism points on
X_0(l) (qf_x0_endos).

Two DIFFERENT matrix encodings of a form coexist here -- don't mix them:
qf_to_mat is the Gram matrix ((2a, b), (b, 2c)) of the bilinear form (used by
the SL2 action act_qf); qf_2_mat is the matrix of multiplication by a*tau on
the basis (1, tau) of the lattice (inverse: mat_2_qf).
"""

from functools import lru_cache

from nt import discfac, gcd
from alg_classes import MatrixElement, Mat_n_Z

M2Z = Mat_n_Z(2)

                                    ###################
                                    # Quadratic Forms #
                                    ###################

### Basics

def qf_ev(qf,xy):
    a,b,c = qf
    x,y = xy
    return a*x*x+b*x*y+c*y*y

def qf_evs_inrange(qf:tuple[int,int,int],m:int):
    a,b,c = qf
    return list({qf_ev(qf,(x,y)) for x in range(-m,m+1) for y in range(-m,m+1)})


def qf_in_fundom(qf:tuple[int,int,int])->bool:
    a,b,c = qf
    if abs(b)< a and a < c:
        return True
    elif b == a and a < c:
        return True
    elif c == a and b >= 0 and b <= a:
        return True
    return False

def qf_gcd(qf:tuple[int,int,int])->int:
    a,b,c = qf
    return gcd(a,gcd(b,c))

def qf_is_prim(qf:tuple[int,int,int])->bool:
    return qf_gcd(qf)==1

def qf_make_prim(qf:tuple[int,int,int])->tuple[int,int,int]:
    a,b,c = qf
    g = gcd(a,gcd(b,c))
    if g > 1:
        a,b,c = a//g,b//g,c//g
    return (a,b,c)

def qf_disc(qf:tuple[int,int,int])->tuple[int,int,int]:
    a,b,c = qf_make_prim(qf)
    return b*b-4*a*c


########################
# Modular group action #
########################

def qf_to_mat(qf:tuple[int,int,int])->MatrixElement:
    a,b,c = qf
    return MatrixElement(((2*a,b),(b,2*c)), M2Z)

def mat_to_qf(m:MatrixElement)->tuple[int,int,int]:
    if not isinstance(m, MatrixElement) or m.n != 2:
        raise TypeError('Input should be 2x2 integer matrix')
    arr = m.vec
    a,b1,b2,c = arr[0][0],arr[0][1],arr[1][0],arr[1][1]
    if b1 != b2:
        raise ValueError('Input should be symmetric matrix')
    if a % 2 != 0 or c % 2 != 0:
        raise ValueError('Diagonal entries should be even')
    return (a//2, b1, c//2)

def act_qf(qf:tuple[int,int,int], m:MatrixElement):
    qfm = qf_to_mat(qf)
    madj = m.adjugate
    qfm_new = madj.transpose * qfm * madj
    return mat_to_qf(qfm_new)


def _qf_reduce(qf:tuple[int,int,int])->tuple[int,int,int]:
    """SL2(Z) reduction of a form into the fundamental domain, with the
    transformed form written out directly -- no MatrixElement arithmetic.  (An
    older matrix-tracking twin, qf_to_fun_dom, validated every 2x2 product
    through the generic ring classes -- millions of checks per rigid-l-set
    search -- and was removed in the 2026-08 cleanup; see tag pre-cleanup.)"""
    while not qf_in_fundom(qf):
        a,b,c = qf
        if a > c:
            qf = (c, -b, a)                              # S
        elif a < abs(b):
            k = b // (2*a)
            if b % (2*a) >= a:
                k += 1
            qf = (a, b - 2*a*k, a*k*k - b*k + c)         # T^-k
        elif a + b == 0:
            qf = (a, b + 2*a, a + b + c)                 # T
        elif a == c and b < 0:
            qf = (c, -b, a)                              # S
        else:
            return qf
    return qf


@lru_cache(maxsize=1<<18)
def qf_mod_gamma(qf:tuple[int,int,int])->tuple[int,int,int]:
    return qf_make_prim(_qf_reduce(qf))

#######################################
# Generating lists of quadratic forms #
#######################################

def get_qfs_all(d:int):
    reps_found = []
    if d % 4 > 1 or d >= 0:
        return reps_found
    b = d % 4
    while 3*b*b <= abs(d):
        num = (b*b-d)//4
        a = b
        while a*a <= num:
            if a == 0:
                a += 1
            if num % a == 0:
                c = num // a
                if qf_in_fundom((a,b,c)):
                    reps_found.append(qf_make_prim((a,b,c)))
                if b != 0 and qf_in_fundom((a,-b,c)):
                    reps_found.append(qf_make_prim((a,-b,c)))
            a += 1
        b += 2
    return reps_found

def get_qfs_strict(d:int):
    return [qf for qf in get_qfs_all(d) if qf_disc(qf)==d]

def qfs_ordered_by_cond(d):
    qfs = get_qfs_all(d)
    qfs.sort(key = lambda qf:discfac(qf_disc(qf))[1])
    return qfs

def class_group_id(d:int):
    if d % 4 > 1:
        raise ValueError(f'{d} is not a discriminant')
    else:
        return (1,d%4,-(d//4))

def class_group_inv(qf:tuple[int,int,int])->tuple[int,int,int]:
    a,b,c = qf_mod_gamma(qf)
    return qf_mod_gamma((a,-b,c))

#############
# Isogenies #
#############
def fricke_inv(qf:tuple[int,int,int],l:int)->tuple[int,int,int]:
    a,b,c = qf
    return (l*l * c, -l*b, a)

def gamma_0_coset_reps(p:int)->list[MatrixElement]:
    return [MatrixElement(((1,0),(0,1)), M2Z)] + [MatrixElement(((0,-1),(1,a)), M2Z)
            for a in range(-(p//2), (p//2)+(p%2))]

def gamma_0_orb(qf:tuple[int,int,int],l:int)->list[tuple[int,int,int]]:
    return [act_qf(qf,m) for m in gamma_0_coset_reps(l)]


def qf_parents(qf:tuple[int,int,int],l:int):
    d = qf_disc(qf)
    return [qf0 for qf0 in _qf_isogs_up_cached(tuple(qf),l) if qf_disc(qf0)>d]

def qf_2_mat(qf):
    a,b,c = qf
    return [[0,-c],[a,-b]]
def mat_2_qf(m):
    a,b,c,d = m[0][0],m[0][1],m[1][0],m[1][1]
    s = c//abs(c)
    return (s*c,s*(a-d),-s*b)


### Computing isogeny codomains
# Two cached enumerations of the index-l sublattice forms:
#   * _qf_isogs_up_cached: only the IMPRIMITIVE (content-l) sublattice forms,
#     i.e. the horizontal + ascending codomains (roots of qf0 mod l).  At most
#     3 entries -- the hot path for the volcano pipeline (cycle walks, rigid
#     searches, parents), so it stays as cheap as it always was.
#   * _qf_isogs_all_cached: ALL l+1 sublattices, with multiplicity, so the
#     descending codomains (disc l^2 d) are included.  qf_isogs exposes this:
#     its output ALWAYS has exactly l+1 elements.
# The content-l entries of the full enumeration reduce to exactly the upward
# ones (reduction commutes with scaling and qf_mod_gamma strips content), so
# the two agree on the horizontal/ascending values.

@lru_cache(maxsize=1<<17)
def _qf_isogs_up_cached(qf0,l):
    a,b,c = qf0
    qfls = []
    if c % l == 0:
        qfls.append(qf_mod_gamma((a*l,b,c//l)))
    for t in range(l):
        qt = a+b*t+c*t*t
        if qt % l == 0:
            at = qt//l
            bt = (b+2*c*t)
            ct = c*l
            qfls.append(qf_mod_gamma((at,bt,ct)))
    return tuple(qfls)

@lru_cache(maxsize=1<<17)
def _qf_isogs_all_cached(qf0,l):
    a,b,c = qf0
    qfls = [qf_mod_gamma((a*l*l,b*l,c))]                    # x -> l x
    for t in range(l):
        qt = a+b*t+c*t*t
        qfls.append(qf_mod_gamma((qt,(b+2*c*t)*l,c*l*l)))   # y -> t x + l y
    return tuple(qfls)

def qf_isogs(qf0,l):
    """All l+1 codomains of the degree-l isogenies from qf0, with multiplicity:
    descending entries have disc l^2 d, horizontal disc d, ascending disc d/l^2."""
    return list(_qf_isogs_all_cached(tuple(qf0),l))

def qf_isogs_hor(qf0,l):
    d = qf_disc(qf0)
    return [qf for qf in _qf_isogs_up_cached(tuple(qf0),l) if qf_disc(qf)==d]


def qf_isogs_des(qf0,l):
    d = qf_disc(qf0)
    return [qf for qf in qf_isogs(qf0,l) if qf_disc(qf)<d]


def qf_isog_parent(qf,l):
    d,c = discfac(qf_disc(qf))
    if c % l != 0:
        raise ValueError(f'No parents in degree {l}')
    else:
        cands = [qf0 for qf0 in _qf_isogs_up_cached(tuple(qf),l)
                 if discfac(qf_disc(qf0))[1]*l == c]
        assert len(cands)==1
        return cands[0]

# Computing isogeny cycles (horizontal walks: the upward enumeration, so the
# descending codomains now in qf_isogs do not derail the cycle)
def qf_isog_cycle(qf0,l):
    cyc = list(_qf_isogs_up_cached(tuple(qf0),l))
    if len(cyc)==1:
        return [qf0,cyc[0]]
    elif len(cyc)>2:
        raise ValueError('Too many isogenies')
    cyc = [qf0,cyc[0]]
    nextbatch = [qf for qf in _qf_isogs_up_cached(cyc[-1],l) if qf not in cyc]
    while len(nextbatch)>0:
        cyc.append(nextbatch[0])
        nextbatch = [qf for qf in _qf_isogs_up_cached(cyc[-1],l) if qf not in cyc]
    return cyc

def qf_isog_cycle_power(qf0,lk):
    l,k = lk
    if k < 0:
        return qf_sibs(qf0,l)
    elif k == 0:
        return [qf0]
    cyc = qf_isog_cycle(qf0,l)
    if k == 1:
        return cyc
    n = len(cyc)
    m = gcd(n,k)
    nm = n//m
    return [cyc[(k*i) % n] for i in range(nm)]

def qf_sibs(qf0:tuple[int,int,int],l:int):
    sibs = qf_isogs_des(qf_isog_parent(qf0,l),l)
    return [qf0]+[qf for qf in sibs if qf != qf0]


##########
# X_0(l) #
##########


def qf_x0_endos(qf:tuple[int,int,int], l:int)->list[tuple[int,int,int]]:
    qf0 = qf_mod_gamma(qf)
    return [qf1 for qf1 in gamma_0_orb(qf0,l) if qf_mod_gamma(fricke_inv(qf1,l))==qf0]


