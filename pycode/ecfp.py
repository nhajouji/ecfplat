from nt import quad_rec,discfac,divisors
from ringclasses import IntegerSquareMatrix
from modularpolynomials import *
from qfs import *
import json
from pathlib import Path

_DATA_DIR = Path(__file__).parent / 'data'

with open(_DATA_DIR / 'ssfp_pc_bij.json', 'r') as f:
    ssfpdata_raw = json.load(f)
ss_precomputed_dictionary = {int(p):ssfpdata_raw[p]for p in ssfpdata_raw}

def get_precomputed_ssdict(p):
    sspcext = {5:{(0,1):(1,0,5),(0,-1):(2,2,3)},7:{(6,1):(1,0,7),(6,-1):(1,1,2)},
           11:{(1,-1):(1,1,3),(1,1):(1,0,11),(0,1):(3,2,4),(0,-1):(3,-2,4)},
           13:{(5,1):(1,0,13),(5,-1):(2,2,7)},
           19:{(18,-1):(1,1,5),(18,1):(1,0,19),(7,1):(4,2,5),(7,-1):(4,-2,5)}}
    if p in sspcext:
        return sspcext[p]
    if p not in ss_precomputed_dictionary:
        raise ValueError(f'No data for p={p}')
    data = ss_precomputed_dictionary[p]
    dic = {}
    for pair in data:
        js = (pair[0][0],pair[0][1])
        qf = (pair[1][0],pair[1][1],pair[1][2])
        dic[js]=qf
    return dic

def sig_to_qfs_dict(a:int,p:int,j0=None):
    if a == 0:
        if p in ss_precomputed_dictionary:
            return get_precomputed_ssdict(p)
        else:
            raise ValueError('No precomputed data for this prime')
    d = a*a-4*p
    if d in hilb_polys_dict and len(hilb_polys_dict[d])==2:
        qf = get_qfs_strict(d)[0]
        j0 = trfr_to_js(a,p)[0]
        return {j0:qf}
    if d >= 0:
        raise ValueError('a is not a trace of Frobenius of ordinary curve mod p')
    leaves = trfr_to_leaves(a,p)
    if j0 == None:
        j0 = leaves[0]
    if j0 not in leaves:
        raise ValueError('j invariant provided does not work')
    qfdata = d_to_ssl_cycle_data(d)
    ml = max([len(qfdata[qf]) for qf in qfdata])
    qf0 = [qf for qf in qfdata if len(qfdata[qf])==ml][0]
    l0 = [l for l in qfdata[qf0] if isinstance(qfdata[qf0][l],list)][0]
    qf_cycle = qfdata[qf0][l0]
    j_cycle = all_ssl_cycles_from_jp(j0,p)[l0]
    qf_to_j = {}
    j_to_qf = {}
    for qj in zip(qf_cycle,j_cycle):
        qf,j = qj
        qf_to_j[qf]= j
        j_to_qf[j]=qf
    for qf1 in qfdata:
        j1 = qf_to_j[qf1]
        for l in qfdata[qf1]:
            if l !=l0:
                qf2 = qfdata[qf1][l]
                j2s = list({j2 for j2 in fp_isog_codomains(j1,l,p) if j2 not in j_to_qf})
                if len(j2s)==1:
                    j2 = j2s[0]
                    qf_to_j[qf2]=j2
                    j_to_qf[j2]=qf2
    c= discfac(a**2-4*p)[1]
    if c == 1:
        return j_to_qf
    for l in atkin_polys_dict:
        while c % l == 0:
            newassignments = {}
            for j0 in j_to_qf:
                j1s = list({j1 for j1 in fp_isog_codomains(j0,l,p) if 
                            j1 not in j_to_qf and j1 not in newassignments})
                if len(j1s)==1:
                    j1 = j1s[0]
                    q0 = j_to_qf[j0]
                    q1s = qf_parents(q0,l)
                    if len(q1s)==1:
                        q1 = q1s[0]
                        newassignments[j1]=q1
            j_to_qf.update(newassignments)
            c = c//l
    return j_to_qf


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

## This is a helper function for the next two
def cubic_qrs(fg:tuple[int,int],p:int)->list[int]:
    f,g = fg
    return [quad_rec((x**3+f*x+g)%p,p) for x in range(p)]

# This returns trace of frobenius, together with the number of points of order 2
# The extra data is needed for supersingular curves
def trace_and_2tor(fg:tuple[int,int],p:int)->tuple[int,int]:
    qrs = cubic_qrs(fg,p)
    return -sum(qrs),len([s for s in qrs if s == 0])

def trace_frob(fg:tuple[int,int],p:int)->int:
    f,g = fg
    return - sum([quad_rec(x**3+f*x+g,p) for x in range(p)])

# Takes a prime as input and returns list of discriminants
# of all endomorphism rings mod p
def supp(p:int)->list[int]:
    ds = []
    a = 0
    while a*a - 4*p<0:
        ds.append(a*a-4*p)
        a+=1
    return ds

def support_closed(p:int)->list[int]:
    suppbd = supp(p)
    supp_clos = []
    for d in suppbd:
        d0,c = discfac(d)
        supp_clos+=[d0*(c0**2) for c0 in divisors(c)]
    ds = list(set(supp_clos))
    ds.sort(key=abs)
    return ds

def supp_in_hilbdb(p:int)->dict:
    ds = support_closed(p)
    return {d:hilb_polys_dict[d] for d in hilb_polys_dict if d in ds}

def endo_db_check(j:int,p:int)->list[int]:
    cands = supp_in_hilbdb(p)
    return [d for d in cands if poly_eval_mod(cands[d][::-1],j,p)==0]


def fp_isog_codomains(j:int,l:int,p:int):
    y0s = atk_at_j(j,l,p).mod(p).find_roots_BrFo()
    atkin_linear = atk_poly_a(l,p).mod(p)
    return [(atkin_linear.eval(y0)-j)%p for y0 in y0s]

##
# This returns cycles of j invariants
# The primes l are chosen so that they do not divide the conductor - 
# this means any two js in the same cycle have the same endomorphism ring
def all_ssl_cycles_from_jp(j0:int,p:int):
    if (j0 == 0 and p % 3 == 1) or ((j0-1728)% p == 0 and p % 4 == 1):
        return {l:[j0] for l in atkin_polys_dict} 
    fg0 = j_to_fg(j0,p)
    a = trace_frob(fg0,p)
    d = a*a-4*p
    d0, c = discfac(d)
    cycles = {}
    for l in atkin_polys_dict:
        if c % l != 0 and quad_rec(d0,l)>=0:
            cycle = []
            nextbatch = [j0]
            while len(nextbatch)>0:
                j1 = nextbatch[0]
                cycle.append(j1)
                nextbatch = [j2 for j2 in fp_isog_codomains(j1,l,p) if j2 not in cycle]
            cycles[l] = cycle
    return cycles
    

#This returns all j-invariants related to j0 by an isogeny of degree l
#It summarizes output of last one, and is needed to deal with noncyclic class groups
#If we don't have a generator, this will contain the codomains of isogenies we found
#and will be bigger than any single orbit in the noncyclic case we hope
def j0_to_j1s_in_sslcycs(j0:int,p:int):
    dic = all_ssl_cycles_from_jp(j0,p)
    codoms = []
    for l in dic:
        codoms+=dic[l]
    j1s = list(set(codoms))
    return j1s
    
def trfr_to_js(a:int,p:int):
    js= [j for j in range(1,p) if 
            (j-1728)%p!=0 and abs(trace_frob(j_to_fg(j),p))==abs(a)]
    d = discfac(a*a-4*p)[0]
    if d==-3 or (d%p == 0 and p % 3 == 2):
        js.append(0)
    if d == -4 or (d % p == 0 and p % 4 == 3):
        js.append(1728%p)
    return js

def trfr_to_models(a:int,p:int):
    d,c = discfac(a**2-4*p)
    j_to_fg_dict = {}
    models = []
    # We need a nonsquare mod p 
    ns = p-1
    while ns>1 and quad_rec(ns,p)>0:
        ns-=1
    if d == -3:
        fg0 = (0,1)
        tr0 = trace_frob(fg0,p)
        if tr0 == a:
            j_to_fg_dict[0]=fg0
        elif tr0 == -a:
            j_to_fg_dict[0] = (0,pow(ns,3,p))
        else:
            gn = p-2
            while pow(gn,p//3,p)==1:
                gn-=1
            tr1 = trace_frob((0,gn),p)
            if tr1 == a:
                j_to_fg_dict[0] = (0,gn)
            elif tr1 == -a:
                j_to_fg_dict[0] = (0,(pow(ns,3,p)*gn)%p)
            else:
                tr2 = trace_frob((0,pow(gn,2,p)),p)
                if tr2 == a:
                    j_to_fg_dict[0]=(0,pow(gn,2,p))
                else:
                    j_to_fg_dict[0] = (0,(pow(ns,3,p)*pow(gn,2,p))%p)
    elif d == -4:
        fg0 = (1,0)
        tr0 = trace_frob(fg0,p)
        if tr0 == a:
            j_to_fg_dict[1728%p]=fg0
        elif tr0 == -a:
            j_to_fg_dict[1728%p]=(pow(ns,2,p),p)
        else:
            tr1 = trace_frob((ns,0),p)
            if tr1 == a:
                j_to_fg_dict[1728%p] = (ns,0)
            else:
                j_to_fg_dict[1728%p] = (pow(ns,3,p),0)
    elif a %p == 0:
        if p % 3 == 2:
            models+=[(0,1),(0,ns)]
        if p % 4 == 3:
            models+=[(p-1,0),(1,0)]
    for j0 in range(1,p):
        if (j0-1728)%p!=0:
            fg0 = j_to_fg(j0,p)
            tr0 = trace_frob(fg0,p)
            if tr0 == a:
                if a == 0:
                    f0, g0 = fg0
                    j_to_fg_dict = [fg0,((f0*(ns**2))%p,(g0*(ns**3))%p)]
                else:
                    j_to_fg_dict[j0] = fg0
                models.append(fg0)
            elif tr0 == -a:
                f0,g0 = fg0
                j_to_fg_dict[j0]=((f0*(ns**2))%p,(g0*(ns**3))%p)
    return j_to_fg_dict


def trfr_to_aec_data(a:int,p:int):
    data = {j:{} for j in trfr_to_js(a,p)}
    d0,c = discfac(a*a-4*p)
    if c == 1:
        return data
    for l in atkin_polys_dict:
        if c % l == 0:
            for j in data:
                data[j][l]= fp_isog_codomains(j,l,p)
    return data
    

def check_leaf(jdata:dict)->bool:
    ls = [l for l in jdata]
    if len(ls) == 0 or len(jdata[ls[0]])==0:
        return True
    for l in ls:
        if len(jdata[l])!= 1:
            return False
    return True

def trfr_to_leaves(a:int,p:int):
    data = trfr_to_aec_data(a,p)
    return [j for j in data if check_leaf(data[j])]
    

def get_j_to_qfs_dict(a:int,p:int,j0=None):
    d = a*a-4*p
    if d in hilb_polys_dict and len(hilb_polys_dict[d])==2:
        qf = get_qfs_strict(d)[0]
        j0 = trfr_to_js(a,p)[0]
        return {j0:qf}
    if d >= 0 or a == 0:
        raise ValueError('a is not a trace of Frobenius of ordinary curve mod p')
    leaves = trfr_to_leaves(a,p)
    if j0 == None:
        j0 = leaves[0]
    if j0 not in leaves:
        raise ValueError('j invariant provided does not work')
    qfdata = d_to_ssl_cycle_data(d)
    ml = max([len(qfdata[qf]) for qf in qfdata])
    qf0 = [qf for qf in qfdata if len(qfdata[qf])==ml][0]
    l0 = [l for l in qfdata[qf0] if isinstance(qfdata[qf0][l],list)][0]
    qf_cycle = qfdata[qf0][l0]
    j_cycle = all_ssl_cycles_from_jp(j0,p)[l0]
    qf_to_j = {}
    j_to_qf = {}
    for qj in zip(qf_cycle,j_cycle):
        qf,j = qj
        qf_to_j[qf]= j
        j_to_qf[j]=qf
    for qf1 in qfdata:
        j1 = qf_to_j[qf1]
        for l in qfdata[qf1]:
            if l !=l0:
                qf2 = qfdata[qf1][l]
                j2s = list({j2 for j2 in fp_isog_codomains(j1,l,p) if j2 not in j_to_qf})
                if len(j2s)==1:
                    j2 = j2s[0]
                    qf_to_j[qf2]=j2
                    j_to_qf[j2]=qf2
    c= discfac(a**2-4*p)[1]
    if c == 1:
        return j_to_qf
    for l in atkin_polys_dict:
        while c % l == 0:
            newassignments = {}
            for j0 in j_to_qf:
                j1s = list({j1 for j1 in fp_isog_codomains(j0,l,p) if 
                            j1 not in j_to_qf and j1 not in newassignments})
                if len(j1s)==1:
                    j1 = j1s[0]
                    q0 = j_to_qf[j0]
                    q1s = qf_parents(q0,l)
                    if len(q1s)==1:
                        q1 = q1s[0]
                        newassignments[j1]=q1
            j_to_qf.update(newassignments)
            c = c//l
    return j_to_qf


## Frobenius matrix
# The following gives a matrix that represents
# the action of Frobenius on the lattice,
# relative to the basis 1, tau.
def frobmat(ap:tuple[int],abc:tuple[int])->IntegerSquareMatrix:
    t,p = ap
    a,b,c = qf_make_prim(abc)
    df,cf = discfac(t*t-4*p)
    dt,ct = discfac(b*b-4*a*c)
    if dt!= df or cf % ct != 0:
        raise ValueError('incompatible discriminants') 
    cft = cf // ct
    trdiff = t+b*cft
    if trdiff%2 != 0:
        raise ValueError('Check trace')
    t0 = trdiff//2
    return IntegerSquareMatrix([[t0,-a*c*cft],[cft,t0-b*cft]])

def kernel_gen_cyc(mat:list[list[int]])->dict:
    m00 = mat[0][0]
    m01 = mat[0][1]
    m10 = mat[1][0]
    m11 = mat[1][1]
    l = m00*m11 - m01*m10
    l0 = gcd(gcd(m00,l),gcd(m01,l))
    l1 = gcd(gcd(m10,l),gcd(m11,l))
    if gcd(l0,l1)>1:
        return 'Check gcds'
    n0 = hall_multiplier(l0,l)
    # Multiplying u by n0 gives us a point that generates a subgroup
    # of index l0c; the order of the point is l/l0c
    l0c = n0*l0
    # We need a point of order l0c, and index l/l0c
    # We already have a point of order l/l1
    n1 = l//(l0c*l1)
    # v already has order l/l1
    # (l/l1) *v = 0, so (l/l1)/(l0c) will have order l0c
    xg = -(n0*m01+n1*m11)%l
    yg = (n0*m00+n1*m10)%l
    # Note that gcd(xg,yg, l) should now be equal to 1.
    # If gcd(xg,yg) is not 1, we can factor it out
    gxy = gcd(xg,yg)
    if gcd(gxy,l)!= 1:
        return 'Something went wrong'
    if gxy == 0:
        return 'Something wrong'
    return {(xg//gxy,yg//gxy):l}


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

# MW generators

def mw_gens(ap:tuple[int],abc:tuple[int],n:int)->dict:
    fmat = frobmat(ap,abc)
    cmat,m = (fmat**n - (fmat**0)).gcdfac()
    mat = IntegerSquareMatrix(cmat.mat).mat
    cdet = mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1]
    if abs(cdet) == 1:
        return {(1,0):m,(0,1):m}
    elif m == 1:
        return kernel_gen_cyc(cmat.mat)
    else:
        return divide_cyclic_gen(kernel_gen_cyc(cmat.mat),m)

    

######
# Old #
# ###    
def get_endo_disc_cands(fg:tuple[int,int],p)->list[int]:
    j0 = fg_to_j(fg,p)
    hds = supp_in_hilbdb(p)
    db_check = [d for d in hds if
                poly_eval_mod(hds[d][::-1],j0,p)==0]
    if len(db_check)>0:
        return db_check
    f,g = fg
    if (f*g) %p != 0:
        a, t2 = trace_and_2tor(fg,p)
        d = a*a-4*p
        d0, c = discfac(d)
        if c == 1:
            return [d]
        elif c == 2:
            if t2==1:
                return [d]
            else:
                return [d0]
        cands = [d0*(c0**2) for c0 in divisors(c)]
        # If the conductor of Z[phi] is even, we can use 2torsion
        # to gain information about the order of 2 in the conductor of O_E
        if c % 2 == 0:
            cands1 = []
            cands2 = []
            for d0 in cands:
                if (d//d0)%2 == 1:
                    cands1.append(d0)
                else:
                    cands2.append(d0)
            if t2 == 1:
                cands = cands1
            else:
                cands= cands2
        #We've already ruled out some discriminants -
        #if they happen to appear in the set of candidates
        #we can discard them immediately
        cands = [d1 for d1 in cands if d1 not in hds]
        return cands
    elif f == 0:
        if p % 3 == 1:
            return [-3]
        else:
            return [-4*p]
    else:
        if p % 4 == 1:
            return [-4]
        qr = quad_rec(-f,p)
        if qr == 1:
            return [p]
        else:
            return [4*p]
    

def j_to_disc(j:int,p:int)->int:
    if j == 0:
        if p % 3 == 1:
            return -3
        else:
            return -4*p
    elif (j-1728)%p == 0:
        if p % 4 == 1:
            return - 4
        else:
            return -4*p
    else:
        f = -3*j*(j-1728)
        g = 2*j*(j-1728)**2
        a = trace_frob((f,g),p)
        return a**2 - 4*p
    