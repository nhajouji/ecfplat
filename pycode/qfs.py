from nt import discfac, gcd
from identities import *
from ringclasses import IntegerSquareMatrix
from modularpolynomials import *

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

def qf_to_mat(qf:tuple[int,int,int])->IntegerSquareMatrix:
    a,b,c = qf
    return IntegerSquareMatrix([[2*a,b],[b,2*c]])

def mat_to_qf(m:IntegerSquareMatrix)->tuple[int,int,int]:
    if not isinstance(m,IntegerSquareMatrix) or m.dim!=2:
        raise TypeError('Input should be 2x2 integer matrix')
    arr = m.mat
    a,b1,b2,c = arr[0][0],arr[0][1],arr[1][0],arr[1][1]
    if b1!= b2:
        raise ValueError('Input should be symmetric matrix')
    if a % 2 !=0 or c %2 != 0:
        raise ValueError('Diagonal entries should be even')
    return (a//2,b1,c//2)

def act_qf(qf:tuple[int,int,int],m:IntegerSquareMatrix):
    qfm = qf_to_mat(qf)
    tm = m.trace()
    minv = IntegerSquareMatrix([[tm,0],[0,tm]])+(-1)*m
    qfm_new = minv.transpose()*qfm*minv
    return mat_to_qf(qfm_new)


def qf_to_fun_dom(qf:tuple)->tuple:
    d = qf_disc(qf)
    if d >=0:
        raise ValueError('Discriminant must be negative')
    matrix = IntegerSquareMatrix([[1,0],[0,1]])
    while not qf_in_fundom(qf):
        a,b,c = qf
        if a>c:
            m0 = IntegerSquareMatrix([[0,-1],[1,0]])
            matrix = m0*matrix
            qf = act_qf(qf,m0)
        elif a < abs(b):
            k = b//(2*a)
            if b % (2*a) >= a:
                k+=1
            m0 = IntegerSquareMatrix([[1,k],[0,1]])
            matrix =  m0*matrix
            qf = act_qf(qf,m0)
        elif a+b==0:
            m0 = IntegerSquareMatrix([[1,-1],[0,1]])
            matrix =  m0*matrix
            qf = act_qf(qf,m0)
        elif a==c and b<0:
            m0 = IntegerSquareMatrix([[0,-1],[1,0]])
            matrix =  m0*matrix
            qf = act_qf(qf,m0)
        else:
            return qf,matrix
    return qf, matrix

def qf_mod_gamma(qf:tuple[int,int,int])->tuple[int,int,int]:
    return qf_make_prim(qf_to_fun_dom(qf)[0])

#######################################
# Generating lists of quadratic forms #
#######################################

def get_qfs_all(d:int):
    reps_found = []
    # First we check that d is indeed a discriminant;
    # if it isn't, we simply return an empty list because there are no associated lattices
    if d % 4 > 1 or d >= 0:
        return reps_found
    b = d % 4
    while 3*b *b <= abs(d):
        num = (b*b-d)//4
        a = b
        while a *a <= num:
            if a == 0:
                a+=1
            if num % a == 0:
                c = num//a
                if qf_in_fundom((a,b,c)):
                    reps_found.append(qf_make_prim((a,b,c)))
                if b!= 0 and qf_in_fundom((a,-b,c)):
                    reps_found.append(qf_make_prim((a,-b,c)))
            a+=1
        b+=2
    return reps_found

def get_qfs_strict(d:int):
    return [qf for qf in get_qfs_all(d) if qf_disc(qf)==d]


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
    return (l*l * c, -l*b,a)

def gamma_0_coset_reps(p:int)->list[IntegerSquareMatrix]:
    return [IntegerSquareMatrix([[1,0],[0,1]])]+[IntegerSquareMatrix([[0,-1],[1,a]])
           for a in range(-(p//2),(p//2)+(p%2))]

def gamma_0_orb(qf:tuple[int,int,int],l:int)->list[tuple[int,int,int]]:
    return [act_qf(qf,m) for m in gamma_0_coset_reps(l)]


def qf_isogenies_all(qf:tuple[int,int,int],l:int,normalize=True):
    isoqfs = [fricke_inv(qf0,l) for qf0 in gamma_0_orb(qf,l)]
    if not normalize:
        return isoqfs
    return [qf_mod_gamma(qf0) for qf0 in isoqfs]

def qf_isogenies_hor(qf:tuple[int,int,int],l:int):
    d = qf_disc(qf)
    return [qf0 for qf0 in qf_isogenies_all(qf,l) if qf_disc(qf0)==d]

def qf_isogenies_down(qf:tuple[int,int,int],l:int):
    d = qf_disc(qf)
    return [qf0 for qf0 in qf_isogenies_all(qf,l) if qf_disc(qf0)<d]

def qf_parents(qf:tuple[int,int,int],l:int):
    d = qf_disc(qf)
    return [qf0 for qf0 in qf_isogenies_all(qf,l) if qf_disc(qf0)>d]

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
        p*=x
    return p

### Computing isogeny codomains
def qf_isogs(qf0,l):
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
    return qfls

def qf_isogs_hor(qf0,l):
    return [qf for qf in qf_isogs(qf0,l) if qf_disc(qf)==qf_disc(qf0)]

def qfs_isogs_int(qfl1,qfl2):
    qf1,l1 = qfl1
    qf2,l2 = qfl2
    qf3s = {qf for qf in qf_isogs_hor(qf1,l1) if qf2 in qf_isogs_hor(qf2,l2)}
    if len(qf3s)==1:
        return list(qf3s)[0]
    else:
        # Something went wrong, either there is no tuple satisfying those conditions
        # or there are multiple tuples - second would be due to an edge case, first should be impossible
        # if pairs are chosen appropriately
        return qf3s
        
# Computing isogeny cycles
def qf_isog_cycle(qf0,l):
    cyc = qf_isogs(qf0,l)
    if len(cyc)==1:
        return [qf0,cyc[0]]
    elif len(cyc)>2:
        raise ValueError('Too many isogenies')
    cyc = [qf0,cyc[0]]
    nextbatch = [qf for qf in qf_isogs(cyc[-1],l) if qf not in cyc]
    while len(nextbatch)>0:
        cyc.append(nextbatch[0])
        nextbatch = [qf for qf in qf_isogs(cyc[-1],l) if qf not in cyc]
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
    sibs = qf_isogenies_down(qf_parents(qf0,l)[0],l)
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
def minv(m:IntegerSquareMatrix)->IntegerSquareMatrix:
    return -m+m.trace()

def find_rrep_g0(m:IntegerSquareMatrix,l:int)->IntegerSquareMatrix:
    reps = gamma_0_coset_reps(l)
    cands= [m0 for m0 in reps if ((m*minv(m0)).mat)[1][0]%l==0]
    if len(cands)!=1:
        raise ValueError('No unique rep')
    else:
        return cands[0]

def qf_to_gamma_0_fd(qf:tuple[int,int,int],l:int)->tuple[tuple[int,int,int],IntegerSquareMatrix]:
    qf0,m = qf_to_fun_dom(qf)
    if m.mat[1][0]%l == 0:
        return qf0, m
    ml = minv(find_rrep_g0(m,l))
    return act_qf(qf,ml),m*ml

def qf_mod_gamma_0(qf:tuple[int,int,int],l:int)->tuple[int,int,int]:
    return qf_to_gamma_0_fd(qf,l)[0]

def qf_x0_endos(qf:tuple[int,int,int],l:int)->list[tuple[int,int,int]]:
    qf0 = qf_mod_gamma(qf)
    return [qf1 for qf1 in gamma_0_orb(qf0,l) if qf_mod_gamma(fricke_inv(qf1,l))==qf0]

def x0_endos_all(p:int)->dict:
    endos_by_trace = {}
    a = 0
    while a*a < 4*p:
        d = a*a-4*p
        qfs = get_qfs_all(d)
        endos_by_trace[a]={qf:qf_x0_endos(qf,p) for qf in qfs}
        a+=1
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
                iso_taus[qf1]=fricke_inv(qf1,l)
    return iso_taus


