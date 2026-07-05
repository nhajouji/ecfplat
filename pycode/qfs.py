from functools import lru_cache

from nt import discfac, gcd
from identities import *
from alg_classes import MatrixElement, Mat_n_Z
# from modularpolynomials import *

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

def qf_to_dc(qf:tuple[int,int,int])->tuple[int,int]:
    return discfac(qf_disc(qf))


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


def qf_to_fun_dom(qf:tuple)->tuple:
    d = qf_disc(qf)
    if d >= 0:
        raise ValueError('Discriminant must be negative')
    matrix = MatrixElement(((1,0),(0,1)), M2Z)
    while not qf_in_fundom(qf):
        a,b,c = qf
        if a > c:
            m0 = MatrixElement(((0,-1),(1,0)), M2Z)
            matrix = m0 * matrix
            qf = act_qf(qf, m0)
        elif a < abs(b):
            k = b // (2*a)
            if b % (2*a) >= a:
                k += 1
            m0 = MatrixElement(((1,k),(0,1)), M2Z)
            matrix = m0 * matrix
            qf = act_qf(qf, m0)
        elif a + b == 0:
            m0 = MatrixElement(((1,-1),(0,1)), M2Z)
            matrix = m0 * matrix
            qf = act_qf(qf, m0)
        elif a == c and b < 0:
            m0 = MatrixElement(((0,-1),(1,0)), M2Z)
            matrix = m0 * matrix
            qf = act_qf(qf, m0)
        else:
            return qf, matrix
    return qf, matrix

def _qf_reduce(qf:tuple[int,int,int])->tuple[int,int,int]:
    """qf_to_fun_dom without the SL2(Z) matrix bookkeeping: same reduction steps
    with the transformed form written out directly, no MatrixElement arithmetic.
    (The matrix-tracking version validated every 2x2 product through the generic
    ring classes -- millions of checks per rigid-l-set search.)"""
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

def prod_tup(t):
    p = 1
    for x in t:
        p *= x
    return p

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

def qf_isogs_asc(qf0,l):
    d = qf_disc(qf0)
    return [qf for qf in _qf_isogs_up_cached(tuple(qf0),l) if qf_disc(qf)>d]

def qf_isogs_des(qf0,l):
    d = qf_disc(qf0)
    return [qf for qf in qf_isogs(qf0,l) if qf_disc(qf)<d]

def qfs_isogs_int(qfl1,qfl2):
    qf1,l1 = qfl1
    qf2,l2 = qfl2
    qf3s = {qf for qf in qf_isogs_hor(qf1,l1) if qf2 in qf_isogs_hor(qf2,l2)}
    if len(qf3s)==1:
        return list(qf3s)[0]
    else:
        return qf3s

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

def cycs_from_ancestors(qf0):
    d, c = discfac(qf_disc(qf0))
    cycs = {}
    if c % 2 == 0:
        sibs = qf_sibs(qf0,2)
        if len(sibs)>1:
            cycs[2] = sibs
    if c % 3 == 0:
        sibs = qf_sibs(qf0,3)
        if len(sibs)<4:
            cycs[3] = sibs
    return cycs


##########
# X_0(l) #
##########

def minv(m:MatrixElement)->MatrixElement:
    return m.adjugate

def find_rrep_g0(m:MatrixElement, l:int)->MatrixElement:
    reps = gamma_0_coset_reps(l)
    cands = [m0 for m0 in reps if (m * minv(m0)).vec[1][0] % l == 0]
    if len(cands) != 1:
        raise ValueError('No unique rep')
    return cands[0]

def qf_to_gamma_0_fd(qf:tuple[int,int,int], l:int)->tuple[tuple[int,int,int], MatrixElement]:
    qf0, m = qf_to_fun_dom(qf)
    if m.vec[1][0] % l == 0:
        return qf0, m
    ml = minv(find_rrep_g0(m, l))
    return act_qf(qf, ml), m * ml

def qf_mod_gamma_0(qf:tuple[int,int,int], l:int)->tuple[int,int,int]:
    return qf_to_gamma_0_fd(qf, l)[0]

def qf_x0_endos(qf:tuple[int,int,int], l:int)->list[tuple[int,int,int]]:
    qf0 = qf_mod_gamma(qf)
    return [qf1 for qf1 in gamma_0_orb(qf0,l) if qf_mod_gamma(fricke_inv(qf1,l))==qf0]

def x0_endos_all(p:int)->dict:
    endos_by_trace = {}
    a = 0
    while a*a < 4*p:
        d = a*a-4*p
        qfs = get_qfs_all(d)
        endos_by_trace[a] = {qf:qf_x0_endos(qf,p) for qf in qfs}
        a += 1
    return endos_by_trace

def iso_taus_x0_l(qf,l):
    qf_reps = gamma_0_orb(qf_mod_gamma(qf),l)
    return [qf1 for qf1 in qf_reps if qf_disc(qf1)==qf_disc(qf)]

def isos_x0_l_all(d,l):
    qfs = get_qfs_all(d)
    iso_taus = {}
    for qf0 in qfs:
        qf1s = gamma_0_orb(qf0,l)
        for qf1 in qf1s:
            if qf_mod_gamma(fricke_inv(qf1,l)) in qfs:
                iso_taus[qf1] = fricke_inv(qf1,l)
    return iso_taus
