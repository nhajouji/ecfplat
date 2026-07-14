import numpy as np
import pandas as pd
import json
from pathlib import Path

from nt import primesBetween,discfac,find_prim_root,quad_rec,gcd,hall_multiplier,axby
from alg_classes import *
from qfs import *
from graph_tools import nbrdata_to_isomat,cycle_from_neighbor_data

M2Z = Mat_n_Z(2)

_DATA_DIR = Path(__file__).parent / 'data'

##################################
# Loading precomputed bijections #
##################################

def strtup_to_tup(s):
    return tuple([int(s0) for s0 in s[1:-1].split(',')])


with open(_DATA_DIR / 'ecqf_ord_pcbij_4_1024.json', 'r') as f:
    ecqf_ord_pcbij_1K_loaded = json.load(f)
ecqf_ord_1K_pc = {strtup_to_tup(aps):{int(ns):tuple(ecqf_ord_pcbij_1K_loaded[aps][ns] )
                                       for ns in ecqf_ord_pcbij_1K_loaded[aps]}
                                       for aps in ecqf_ord_pcbij_1K_loaded}

with open(_DATA_DIR / 'ecqf_ss_pcbij_velu_4_1024.json', 'r') as f:
    ecqf_ss_pcbij_1K_loaded = json.load(f)
ecqf_ss_1K_pc = {int(ps):{strtup_to_tup(jss):tuple(ecqf_ss_pcbij_1K_loaded[ps][jss] )
                                       for jss in ecqf_ss_pcbij_1K_loaded[ps]} 
                                       for ps in ecqf_ss_pcbij_1K_loaded}


def get_aps_pc():
    return [ap for ap in ecqf_ord_1K_pc]

def get_ssps_pc():
    ps = [p for p in ecqf_ss_1K_pc]
    ps.sort()
    return ps

def ap_in_pc_data(ap:tuple[int,int])->bool:
    if ap[0]==0:
        return ap[1] in get_ssps_pc()
    else:
        a,p = ap
        if (abs(a),p) in get_aps_pc():
            return True
        # beyond the eagerly-loaded 1K store: the lazy sqlite ext store
        import ext_store
        return ext_store.has_pair(abs(a), p)

# Misc. tools
def abc_to_tau(abc:tuple[int,int,int])->np.array:
    tau_np= np.roots(abc)[0]
    return np.array([tau_np.real,tau_np.imag])

def abc_to_tau_str(abc):
    a,b,c = qf_make_prim(abc)
    d,cond = discfac(b*b-4*a*c)
    den = 2*a
    g = gcd_list([den,b,cond])
    if g>1:
        den,b,cond = den//g, b//g, cond//g
    if d % 4 == 0 and b % 2 == 0 and den % 2 == 0:
        d = d//4
        b = b//2
        den = den//2
    num_str = f' √{d}'
    if cond > 1:
        num_str = f"{cond}"+num_str
    if b!=0:
        num_str = f"{-b}+"+num_str
    if den == 1:
        return num_str
    if '+' in num_str:
        return '('+num_str+f')/{den}'
    else:
        return num_str+'/'+str(den)

def ec_eq_str_base(fg:tuple[int,int]):
    """'y^2 = x^3 + f x + g' with 0/±1/negative coefficients handled, no 'mod p'."""
    eqs = 'y^2 = x^3'
    f,g = fg
    if f!= 0:
        if f>0:
            if f == 1:
                eqs+= ' + x'
            else:
                eqs+=f' + {f} x'
        else:
            f = -f
            if f == 1:
                eqs+= ' - x'
            else:
                eqs+=f' - {f} x'
    if g != 0:
        if g>0:
            eqs+=f' + {g}'
        else:
            g = -g
            eqs+=f' - {g}'
    return eqs

def ec_eq_str(fg:tuple[int,int],p:int):
    return ec_eq_str_base(fg)+(f' mod {p}')

def export_points(grp:list,filename:str):
    file = open(filename,'a')
    for i in range((len(grp)//3)+1):
        stri = ''
        for vec in grp[3*i:3*(i+1)]:
            stri += '['+str(vec[0])[:10]+','+str(vec[1])[:10]+'],'
        stri+='\n'
        file.write(stri)
    file.close()

def qf_cond(qf:tuple[int,int,int])->int:
    return discfac(qf_disc(qf))[1]

def m2_tup_to_vec(m:tuple[tuple[int]])->tuple[int,int,int,int]:
    return tuple(list(m[0])+list(m[1]))

def m2_tup_det(m:tuple[tuple[int]])->int:
    a,b,c,d = m2_tup_to_vec(m)
    return a*d-b*c


def pts_from_gendic(gendic:dict)->list[tuple[int,int]]:
    if len(gendic)== 0:
        return [(0,0)]
    gen_tups = [v for v in gendic]
    if len(gen_tups)==2 and gendic[gen_tups[0]]==gendic[gen_tups[1]]:
        l = gendic[gen_tups[0]]
        return [(a,b) for a in range(l) for b in range(l)]
    lm = max(gendic.values())
    A = ZnProduct((lm.lm))
    if len(gen_tups)==1:
        v = gen_tups[0]
        g = AbGrElt(v,A)
        return [(m*g).vec for m in range(lm)]
    else:
        v1,v2 = gen_tups[1],gen_tups[2]
        g1,g2 = AbGrElt(v1,A),AbGrElt(v2,A)
        n1,n2 = gendic[v1],gendic[v2]
        return [(x*g1+y*g2).vec for x in range(n1) for y in range(n2)]




######
# ls #
######

def qf_l_order(qf,mm = 64):
    d = qf_disc(qf)
    if qf == class_group_id(d):
        return 1,1
    m = 8
    evs = qf_evs_inrange(qf,m)
    pevs = [p for p in primesBetween(1,max(evs)) if p in evs]
    while len(pevs)==0 and m<=mm:
        m*=2
        evs = qf_evs_inrange(qf,m)
        pevs = [p for p in primesBetween(1,max(evs)) if p in evs]
    if len(pevs)==0:
        raise ValueError('No primes found, try bigger bound')
    l = min(pevs)
    return l,len(qf_isog_cycle(qf,l))

### ECFP

def trace_frob(fg:tuple[int,int],p:int)->int:
    f,g = fg
    return - sum([quad_rec(x**3+f*x+g,p) for x in range(p)])
def fg_to_j(fg:tuple[int,int],char =0):
    f,g = fg
    if f == 0 or char>0 and f%char ==0:
        return 0
    elif g == 0 or char>0 and g % char ==0:
        return 1728
    else:
        f3 = 4*(f**3)
        jnum = 1728 *f3
        jden = f3+27*(g**2)
        if char == 0:
            if jden == 0:
                raise ZeroDivisionError('Singular curve')
            elif jnum % jden == 0:
                return jnum//jden
            else:
                return jnum/jden
        else:
            jnum = jnum % char
            jden = jden % char
            if jden == 0:
                raise ZeroDivisionError('Singular curve')
            jdeninv = pow(jden,-1,char)
            return (jnum*jdeninv)%char


def j_to_fg(j:int,char = 0):
    if j == 0:
        return (0,1)
    elif j == 1728 or char>0 and (j-1728)%char == 0:
        return (1,0)
    else:
        f = -3*j*(j-1728)
        g = 2*j*((j-1728)**2)
        if char == 0:
            return (f,g)
        else:
            return (f % char, g% char)
        
def twist_fg(fg,t,char = 0):
    f,g = fg
    ft, gt = t**2 * f, t**3*g
    if char == 0:
        return (ft,gt)
    else:
        return (ft%char,gt%char)
    

def ecfp_js_to_model(js,ap,nsp):
    a,p = ap
    if a == 0:
        if not isinstance(js,tuple) or len(js)!=2 or abs(js[1])!=1:
            raise ValueError('Supersingular signature should be pair')
        j,s = js
        if (j-1728)%p == 0:
            if s == 1:
                return (1,0)
            else:
                return (p-1,0)
        elif j == 0:
            if s == 1:
                return (0,1)
            else:
                return (0,p-3)
        f,g = j_to_fg(j,p)
        if quad_rec(g,p)==s:
            return (f,g)
        else:
            return twist_fg((f,g),nsp,p) 
    else:
        if isinstance(js,tuple):
            j = js[0]
        else:
            j = js
        if (j*(j-1728))% p != 0:
            fg = j_to_fg(j,p)
            if trace_frob(fg,p)==a:
                return fg
            else:
                return twist_fg(fg,nsp,p)
        xp = find_prim_root(p)
        if j == 0:
            cands = [(0,pow(xp,k,p)) for k in range(6)]
            fga = [fg for fg in cands if trace_frob(fg,p)==a]
            if len(fga) == 0:
                raise ValueError('No model found')
            return fga[0]
        else:
            cands = [(pow(xp,k,p),0) for k in range(4)]
            fga = [fg for fg in cands if trace_frob(fg,p)==a]
            if len(fga) == 0:
                raise ValueError('No model found')
            return fga[0]



                ###################################
                # Mordell-Weil Group Computations #
                ###################################


# We want to compute fixed points of a power of frobenius
# Just need ability to compute kernel of a general matrix

                        # Kernels #

# Step 1: Factor M as gM0, where g is the gcd of the entries of M 
# and M0 is a matrix with relatively prime entries

def m2_tup_gcdfact(m:tuple[tuple[int]]):
    mv = m2_tup_to_vec(m)
    a,b,c,d = mv
    gcdm = gcd(gcd(a,b),gcd(c,d))
    a,b,c,d = tuple([m0//gcdm for m0 in mv])
    return gcdm, ((a,b),(c,d))

# Step 2: This computes the kernel of a matrix with relatively prime entries

def qf_mat_ker_cyc(mat:tuple[tuple[int]])->dict:
    a,b,c,d = m2_tup_to_vec(mat)
    l = a*d-b*c
    l0 = gcd(gcd(a,b),l)
    l1 = gcd(gcd(c,d),l)
    if gcd(l0,l1)!= 1:
        raise ValueError('Check gcds')
    n0 = hall_multiplier(l0,l)
    # Multiplying u by n0 gives us a point that generates a subgroup
    # of index l0c; the order of the point is l/l0c
    l0c = n0*l0
    # We need a point of order l0c, and index l/l0c
    # We already have a point of order l/l1
    n1 = l//(l0c*l1)
    # v already has order l/l1
    # (l/l1) *v = 0, so (l/l1)/(l0c) will have order l0c
    xg = -(n0*b+n1*d)%l
    yg = (n0*a+n1*c)%l
    # Note that gcd(xg,yg, l) should now be equal to 1.
    # If gcd(xg,yg) is not 1, we can factor it out
    gxy = gcd(xg,yg)
    if gcd(gxy,l)!= 1:
        return 'Something went wrong'
    if gxy == 0:
        return 'Something wrong'
    return {(xg//gxy,yg//gxy):l}

# This combines the previous functions - given a matrix, it returns generators of the kernel as a dictionary.
# The length of the dictionary is the number of generators - 0, 1, or 2 depending on trivial, cyclic, general
### The kernel of M0 is cyclic, we compute that first and then divide by the gcd using:


# Step 3: Extend to a basis of the whole group

def divide_cyclic_gen(gen:dict,m:int)->dict:
    v = [v for v in gen][0]
    l = gen[v]
    x,y = v
    if gcd(x,y)>1:
        g =gcd(x,y)
        if gcd(g,l)>1:
            return'The generator has the wrong order'
        v = x//g, y//g
        x,y = v
    v1 = (l*x)%(l*m),(l*y)%(l*m) 
    r,s = axby(v1)
    w = ((l*s) %(l*m),(l*r)%(l*m))
    return {v:m*l,w:m}

# This combines the previous 3 steps to combine the kernel of a matrix.

def qf_mat_ker_gens(mat:tuple[tuple[int]]):
    if m2_tup_det(mat)==1:
        return {}
    m_gcd,m_cyc = m2_tup_gcdfact(mat)
    if m2_tup_det(m_cyc) == 1:
        return {(1,0):m_gcd,(0,1):m_gcd}
    elif m_gcd == 1:
        return qf_mat_ker_cyc(m_cyc)
    return divide_cyclic_gen(qf_mat_ker_cyc(m_cyc),m_gcd) 

                    # Frobenius Matrix #

# This computes a matrix that represents a generator of the endomorphism ring,
# relative to the basis 1, tau of the associated lattice
# The generator is multiplication-by-(a * tau), where a is the leading coef

def qf_to_ERGM_1T(qf:tuple[int,int,int])->MatrixElement:
    a,b,c = qf
    return MatrixElement(((0,-c),(a,-b)),M2Z)

# This computes the matrix that represents multiplication by
#  the root of x^2-ax + p whose imaginary part has sign equal to s
#  on the lattice relative to the ordered basis 1, tau

def qf_ap_FrMat(qf:tuple[int,int,int],ap:tuple[int,int],s=1)->MatrixElement:
    a,p = ap
    d_qf = qf_disc(qf)
    d_ap = a*a-4*p
    dc_qf = discfac(d_qf)
    dc_ap = discfac(d_ap)
    c_qf = dc_qf[1]
    c_ap = dc_ap[1]
    if dc_qf[0]!=dc_ap[0] or c_ap % c_qf!= 0:
        raise ValueError('Incompatible forms')
    I2m = MatrixElement(((1,0),(0,1)),M2Z)
    atau_mat = qf_to_ERGM_1T(qf)
    tau_scalar = s*(c_ap//c_qf)
    trace_diff = a- (tau_scalar * atau_mat.trace)
    assert trace_diff % 2 == 0
    one_scalar = trace_diff//2
    frmat = one_scalar * I2m + tau_scalar * atau_mat
    if abs(frmat.trace)!= abs(a):
        raise ValueError('Something went wrong in computation')
    if frmat.trace == a:
        return frmat
    else:
        return -frmat


def frk_fxp_mat(frmat:MatrixElement,k:int)->tuple[tuple[int,int],tuple[int,int]]:
    return (frmat**k-frmat**0).vec

# Given a frobenius matrix M and an integer k, this returns generators of the kernel of M^k-1
# In general, this computes generators of the group of fixed points of M^k
def frob_to_mw_gens(frmat:MatrixElement,k:int):
    return qf_mat_ker_gens(frk_fxp_mat(frmat,k))

# The following computes Frobenius matrices for all qfs in a given isogeny class
def ap_FrbMats_1T(ap:tuple[int,int],s=1):
    a,p = ap
    d = a*a-4*p
    d0,c_ap = discfac(d)
    qfs = get_qfs_all(d)
    I2m = MatrixElement(((1,0),(0,1)),M2Z)
    qf_to_Fr_mat = {}
    for qf in qfs:
        atau_mat = qf_to_ERGM_1T(qf)
        c_qf = qf_cond(qf)
        assert c_ap % c_qf == 0
        tau_scalar = s*(c_ap//c_qf)
        trace_diff = a- (tau_scalar * atau_mat.trace)
        assert trace_diff % 2 == 0
        one_scalar = trace_diff//2
        frmat_tup = (one_scalar * I2m + tau_scalar * atau_mat).vec
        frmat= MatrixElement(frmat_tup,M2Z)
        qf_to_Fr_mat[qf] = frmat
    return qf_to_Fr_mat

def mw_arr_from_gens(abc:tuple,gens:dict)->np.array:
    den = max(gens.values())
    one = np.array([1,0])
    tau = abc_to_tau(abc)
    pts = [np.array([0,0])]
    for gen in gens:
        x,y = gen
        pts = [(pt0+m*np.array([x,y]))%den 
               for pt0 in pts for m in range(gens[gen])]
    return np.array([pt[0]*one+pt[1]*tau for pt in pts])/den


                #####################
                # Single Curve Data #
                #####################

def ec_look_up(fg:tuple[int,int],p:int)->dict:
    a = trace_frob(fg,p)
    ap0 = (abs(a),p)
    ap = (a,p)
    f,g = fg
    ecd = (4*pow(f,3,p)+27*pow(g,2,p))%p
    if ecd == 0:
        raise ValueError('This is a singular curve')    
    if g == 0:
        s = quad_rec(f,p)
    else:
        s = quad_rec(g,p)
    j = fg_to_j(fg,p) %p
    df = a*a-4*p
    d,cf = discfac(df)
    data = {'ap':ap,'ec_eq':ec_eq_str(fg,p),'coefs':fg,
            'j':j,'s':s,'frob_tr':a,'frob_disc':df,'frob_cond':cf,'endo_fun_disc':d,
            'has_pcqf':False,'qf':None}
    if a == 0:
        if p in get_ssps_pc():
            data['has_pcqf'] = True
            data['qf'] = ecqf_ss_1K_pc[p][(j,s)]
    else:
        if ap0 in get_aps_pc():
            data['has_pcqf'] = True
            data['qf'] = ecqf_ord_1K_pc[ap0][j%p]
    if not data['has_pcqf']:
        return data
    qf = data['qf']
    fr_s = 1
    if a < 0:
        fr_s = -1
    data['FrobMat']= qf_ap_FrMat(qf,ap,s=fr_s)
    tau_arr = abc_to_tau(qf)
    data['tau_str'] = abc_to_tau_str(qf)
    data['tau_xy'] = tuple([np.round(x,3) for x in tau_arr])
    return data

            ##############
            # ECQF Class #
            ##############

# ---------------------------------------------------------------------------
# QFIsogenyClass and ECQFIsogenyClass moved to ecqf.py (2026-07-13) — ecqf is
# now the assembly point the website imports from.
# ---------------------------------------------------------------------------
