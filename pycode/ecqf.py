import itertools
from qfs import *
from nt import discfac,quad_rec,divisors,primesBetween
from identities import *
from ecfp import fp_isog_codomains, trfr_to_js
from modularpolynomials import *

ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
heeg_dics = {}
for d in heeg_js:
    d0 = discfac(d)[0]
    if d0 == d:
        heeg_dics[d]={heeg_js[d]:class_group_id(d)}
    else:
        heeg_dics[d] = {heeg_js[d0]:class_group_id(d0),heeg_js[d]:class_group_id(d)}

def small_bij_check(d):
    if d in heeg_dics:
        return heeg_dics[d]
    else:
        return {}
#########################
# Load precomputed data #
#########################
import json 

def strtup_to_tup(s):
    return tuple([int(s0) for s0 in s[1:-1].split(',')])

with open('ecc/data/ecqf_ord_pcbij_4to256.json', 'r') as f:
    ecqf_ord_pcbij_4to256_loaded = json.load(f)
ecqf_ord_256_pc = {strtup_to_tup(aps):{int(ns):tuple(ecqf_ord_pcbij_4to256_loaded[aps][ns] )
                                       for ns in ecqf_ord_pcbij_4to256_loaded[aps]} 
                                       for aps in ecqf_ord_pcbij_4to256_loaded}



##############
# Misc tools #
##############
def sortpair(ls,ns):
    assert len(ls)==len(ns)
    pairs = [(ls[i],ns[i]) for i in range(len(ls))]
    pairs.sort(key = lambda ln:ln[1],reverse = True)
    return tuple([ln[0] for ln in pairs]),tuple([ln[1] for ln in pairs])

def ext_tup(t,x):
    return tuple(list(t)+[x])

def invert_dict(dict):
    return {dict[x]:x for x in dict}


# Given two dictionaries with the same keys, returns
# a dictionary whose keys are keys of dict, and values are keys of second dict
def compdiv_dics(dic1:dict,dic2:dict)->dict:
    if len(dic1)!=len(dic2):
        raise ValueError('Dictionaries must have same size')
    common_elts = {d for d in dic1 if d in dic2}
    if len(common_elts)<len(dic1):
        raise ValueError('Dictionaries should have same keys')
    dic1inv = {dic1[x]:x for x in dic1}
    return {y:dic2[dic1inv[y]] for y in dic1inv}

###########################
# Obtaining generating ls #
###########################

def scoredic_to_tups(scores:dict):
    ls = [l for l in scores]
    ls.sort(reverse=True)
    return tuple(ls),tuple([scores[l] for l in ls])




def qf_reps_pm(d):
    return [qf for qf in get_qfs_strict(d) if qf[1]>=0]


def disc_to_ssls(d):
    ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
    qfsfound = {l:[] for l in ssprimes}
    qf_reps = qf_reps_pm(d)
    ls = {}
    qf0 = class_group_id(d)
    for qf in qf_reps:
        qfls = [l for l in qf_evs_inrange(qf,10) if l in ssprimes]
        if len(qfls)>0:
            l0 = min(qfls)
            ls[l0]=len(qf_isog_cycle(qf0,l0))
    return ls

def disc_to_ssl_qfs(d):
    ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
    qfsfound = {l:[] for l in ssprimes}
    qf_reps = qf_reps_pm(d)
    ls = {}
    qf0 = class_group_id(d)
    for qf in qf_reps:
        qfls = [l for l in qf_evs_inrange(qf,10) if l in ssprimes]
        if len(qfls)>0:
            for l in qfls:
                ls[qf] = min(qfls)
    return ls

def gen_qfs_d_ls(d,lset):
    qf0 = class_group_id(d)
    qfs = [qf0]
    for l in lset:
        qfs_new = []
        for qf in qfs:
            qfs_new+=qf_isog_cycle(qf,l)
        qfs = list(set(qfs_new))
    return qfs

def score_tuple(ls,scores):
    prod = 1
    for l in ls:
        prod*=scores[l]
    return prod

def check_lset(d,ls):
    cld = clgr_size_gen(d)
    lgr_size = len(gen_qfs_d_ls(d,ls))
    prd = 1
    nlist = []
    for l in ls:
        nl = len(gen_qfs_d_ls(d,[l]))
        prd*=nl
    return [cld,lgr_size,prd]

def minimize_lset(d,ls):
    ls = [l for l in ls]
    n1 = len(gen_qfs_d_ls(d,ls))
    l0s= [ls[:i]+ls[i+1:] for i in range(len(ls))]
    l0s_g =[l0 for l0 in l0s if len(gen_qfs_d_ls(d,l0))==n1]
    while len(l0s_g)>0:
        l0s_g.sort(key = max)
        ls = l0s_g[0]
        l0s= [ls[:i]+ls[i+1:] for i in range(len(ls))]
        l0s_g =[l0 for l0 in l0s if len(gen_qfs_d_ls(d,l0))==n1]
    return tuple(ls)

def get_spanl2s(d):
    l2s = [ln[0] for ln in disc_to_ssls(d).items() if ln[1]==2]
    mx = min(len(l2s),clgr_2len(d)-1)
    if mx < 2:
        return [(l,) for l in l2s]
    all_ls = {}
    for k in range(2,mx+1):
        subs = []
        for lks in itertools.combinations(l2s,k):
            if len(gen_qfs_d_ls(d,list(lks)))==2**k:
                subs.append(lks)
        all_ls[k]=subs
    ks_nz = [k for k in all_ls if len(all_ls[k])>0]
    k = max(ks_nz)
    return all_ls[k]

def clgr_n2k_ssls(d):
    cld = clgr_size_gen(d)
    if cld == 1:
        return [()]
    ldata = disc_to_ssls(d)
    if len(ldata)==0:
        raise ValueError('No data found')
    ml = max(ldata.values())
    lms = [l for l in ldata if ldata[l]==ml]
    if ml == cld:
        return [(l,) for l in lms]
    ls2s = get_spanl2s(d)
    if ml == 2:
        return ls2s
    ln2s = []
    for lm in lms:
        for ls2 in ls2s:
            ln2 = [lm]+list(ls2)
            if len(gen_qfs_d_ls(d,ln2))==cld:
                return [minimize_lset(d,ln2)]
    return []


def disc_search_ssls_by_size(d,k=None):
    cld = clgr_size_gen(d)
    ssls = disc_to_ssls(d)
    ls_all = [l for l in ssls]
    if k == None:
        k = len(minimize_lset(d,ls_all))
    ls_gen = {}
    for ls in itertools.combinations(ls_all,k):
        if len(gen_qfs_d_ls(d,list(ls))) == cld:
            ls_gen[ls]=score_tuple(ls,ssls)
    if len(ls_gen)== 0:
        return []
    m = min(ls_gen.values())
    return [{l:ssls[l] for l in ls} for ls in ls_gen if ls_gen[ls]==m]

def disc_search_ssl_gen(d):
    n2k_ls = clgr_n2k_ssls(d)
    if len(n2k_ls)>0:
        return n2k_ls
    else:
        scores = disc_search_ssls_by_size(d)
        return disc_search_ssls_by_size(d)

def qf_trios_from_ltrio(d,l123):
    qf0 = class_group_id(d)
    l1,l2,l3 = l123
    qfs_l1 = qf_isogs_hor(qf0,l1)
    qfs_l2 = qf_isogs_hor(qf0,l2)
    qfs_l3 = qf_isogs_hor(qf0,l3)
    qftrios = []
    for qf1 in qfs_l1:
        qf1_l2_qfs = qf_isogs_hor(qf1,l2)
        for qf2 in qfs_l2:
            qfs_l123 = [qf for qf in qf_isogs_hor(qf2,l1) if qf in qf1_l2_qfs and qf in qfs_l3]
            if len(qfs_l123)>0:
                qftrios.append([qf1,qf2,qfs_l123[0]])
    return qftrios

def qf_nn_rig_basis(d,ldata=None):
    if ldata == None:
        ldata = disc_rig_ssl_search(d)
    if len(ldata['ls'])!=2 or ldata['l_sum'] == None:
        raise ValueError('Use diff algorithm')
    l1,l2 = ldata['ls']
    l3 = ldata['l_sum']
    trios = qf_trios_from_ltrio(d,(l1,l2,l3))
    if len(trios)== 0:
        raise ValueError('No basis found')
    qf1,qf2,qf3 = tuple(trios[0])
    return {():class_group_id(d),(l1,):qf1,(l2,):qf2,(l1,l2):qf3,l3:qf3}


def qf_rig_ssl_data(d,lscores):
    ls_big= [l for l in lscores if lscores[l]>2]
    ls,ns = scoredic_to_tups(lscores)
    data = {'ls':ls,'ns':ns,'needs_sum':False,'l_sum':None}
    if len(ls_big)<2:
        return data
    data['needs_sum']=True
    qf0 = class_group_id(d)
    qfs = [qf0]
    for l in ls_big:
        qfs_ext = []
        for qf in qfs:
            qfs_ext+=qf_isogs_hor(qf,l)
        qfs = qfs_ext
    qf_to_l_data = disc_to_ssl_qfs(d)
    for qf in qfs:
        if qf in qf_to_l_data:
            data['l_sum'] = qf_to_l_data[qf]
            if len(ls_big)==2:
                rigbasis = qf_nn_rig_basis(d,data)
                data['nn_frame'] = rigbasis
            return data
    return data



def disc_rig_ssl_search(d):
    cld = clgr_size_gen(d)
    if cld == 1:
        return {'ls':(),'ns':(),'needs_sum':False,'l_sum':None}
    ssls = disc_to_ssls(d)
    ls_n2k = clgr_n2k_ssls(d)
    if len(ls_n2k)>0:
        ls = ls_n2k[0]
        n0 = max({ssls[l] for l in ls})
        ns = tuple([n0]+[2 for _ in range(len(ls)-1)])
        return {'ls':ls,'ns':ns,'needs_sum':False,'l_sum':None}
    ls_all = [l for l in ssls]
    k = len(minimize_lset(d,ls_all))
    ls_gen = {}
    for ls in itertools.combinations(ls_all,k):
        if len(gen_qfs_d_ls(d,list(ls))) == cld:
            ls_gen[ls]=score_tuple(ls,ssls)
    if len(ls_gen)== 0:
        return []
    m = min(ls_gen.values())
    if m != cld:
        return [ls for ls in ls_gen if ls_gen[ls]==m]
    ls_gen_l = [l for l in ls_gen]
    for ls in ls_gen_l[:-1]:
        rigdata = qf_rig_ssl_data(d,{l:ssls[l] for l in ls})
        if (rigdata['l_sum']!=None) or (not rigdata['needs_sum']):
            return rigdata
    return qf_rig_ssl_data(d,{l:ssls[l] for l in ls_gen_l[-1]})



###################################
# Reconstructing the Cayley graph #
###################################

def qf_isopwr_ngbrs(qf0,lk):
    cyc = qf_isog_cycle_power(qf0,lk)
    if len(cyc)<2:
        return []
    elif len(cyc)==2:
        return [cyc[1]]
    else:
        return [cyc[k] for k in [1,-1]]

def qf_isopwr_intrs(qfl0,qfl1):
    qf0,l0 = qfl0
    qf1,l1 = qfl1
    return [qf2 for qf2 in qf_isopwr_ngbrs(qf0,l0) if qf2 in qf_isopwr_ngbrs(qf1,l1)]


def fp_isog_cycle(jp,l):
    j,p = jp
    cyc = [j]
    nbrs = fp_isog_codomains(j,l,p)
    if len(nbrs)==0:
        return cyc
    cyc.append(nbrs[0])
    if len(nbrs)==1:
        return cyc
    nextbatch = [j2 for j2 in fp_isog_codomains(cyc[-1],l,p) if j2 not in cyc]
    while len(nextbatch)>0:
        cyc.append(nextbatch[0])
        nextbatch = [j2 for j2 in fp_isog_codomains(cyc[-1],l,p) if j2 not in cyc]
    return cyc

def fp_isog_cycle_power(jp,lk):
    l,k = lk
    cycl = fp_isog_cycle(jp,l)
    n = len(cycl)
    kord = len(set((k*a)%n for a in range(n)))
    return [cycl[(k*a)%n] for a in range(kord)]

def fp_isopwr_ngbrs(jp,lk):
    cyc = fp_isog_cycle_power(jp,lk)
    if len(cyc)<2:
        return []
    elif len(cyc)==2:
        return [cyc[1]]
    else:
        return [cyc[k] for k in [1,-1]]

def fp_isopwr_intrs(jpl0,jpl1):
    jp0,l0 = jpl0
    jp1,l1 = jpl1
    return [j for j in fp_isopwr_ngbrs(jp0,l0) if j in fp_isopwr_ngbrs(jp1,l1)]



#######################
# End to end bijection #
########################

# Currently only works for n2k groups #

def tree_edges_to_ancestors(nbrdata):
    leaves = [v for v in nbrdata if len(nbrdata[v])==1 and v != 0]
    anc_data = {}
    nextbatch = []
    for v0 in leaves:
        v1 = nbrdata[v0][0]
        anc_data[v0]=v1
        nextbatch.append(v1)
    while len(nextbatch)>0:
        currentbatch = nextbatch
        nextbatch = []
        for v in currentbatch:
            v_ancs = [v1 for v1 in nbrdata[v] if v1 not in anc_data]
            if len(v_ancs)==1:
                anc_data[v]=v_ancs[0]
                nextbatch+=v_ancs
    return anc_data



def js_to_rabs(js,p):
    ab_to_js = {}
    for i,j1 in enumerate(js):
        for j2 in js[i:]:
            ab_to_js[((-j1-j2)%p,(j1*j2)%p)] = (j1,j2)
    return ab_to_js

def get_nbr_data(ap,l):
    a,p = ap
    js = trfr_to_js(a,p)
    rabs = js_to_rabs(js,p)
    nbrdata = {j:[] for j in js}
    for x in range(p):
        evx = eval_atk(x,l,p)
        if evx in rabs:
            j1,j2 = rabs[evx]
            nbrdata[j1].append(j2)
            nbrdata[j2].append(j1)
    return nbrdata


def get_ancestor_data_ord(ap):
    a,p = ap
    d,c = discfac(a*a-4*p)
    js = trfr_to_js(a,p)
    if c == 1:
        return {'ancestor_data':{},'leaves':js,'js_all':js}
    elif c == 2:
        if d == -4:
            return {'ancestor_data':{287496%p:1728%p},'leaves':[287496%p],'js_all':js}
        leaves = [j for j in js if quad_rec(j-1728,p)==-1]
        ancs = {}
        nbrs = get_nbr_data(ap,2)
        for j in leaves:
            ancs[j]=nbrs[j][0]
        return {'ancestor_data':{2:ancs},'leaves':leaves,'js_all':js}
    ls = [l for l in primesBetween(1,c+1) if c%l ==0]
    anc_data = {}
    leaf_cands = [j for j in js if j!= 0 and (j-1728)%p !=0]
    for l in ls:
        nbrs_l = get_nbr_data(ap,l)
        anc_data[l]=tree_edges_to_ancestors(nbrs_l)
        leaf_cands = [j for j in leaf_cands if len(nbrs_l[j])==1]
    return {'ancestor_data':anc_data,'leaves':leaf_cands,'js_all':js}





##########################################
# Leaf Bij. algorithms for special cases #
##########################################

###############
# (n,2,2,2..) #
###############

def ecqf_ord_n2_bij_leaves(ap,vertical_iso_data=None,ldata = None):
    a,p = ap
    if a == 0:
        raise ValueError('Use supersingular algorithm')
    d = a*a-4*p
    if vertical_iso_data == None:
        vertical_iso_data = get_ancestor_data_ord(ap)
    leaves = [j for j in vertical_iso_data['leaves'] if j*(j-1728)%p!=0]
    if len(leaves)<4:
        # Any bijection works 
        # Will return a random one to avoid dealing with weird issues with descendents of j = 0,1728
        qfs = get_qfs_strict(d)
        assert len(qfs)== len(leaves)
        return {ji:qfs[i] for i,ji in enumerate(leaves)}
    j0 = leaves[0]
    # Look up degrees that are guaranteed to give connected graph
    if ldata == None:
        ldata = disc_rig_ssl_search(d)
    if ldata['needs_sum']:
        raise ValueError('Use different algorithm')
    j_to_qf = {j0:class_group_id(d)}
    qf_to_j = {class_group_id(d):j0}
    ls = ldata['ls']
    for lk in ls:
        newassignments = {}
        for qf in qf_to_j:
            j = qf_to_j[qf]
            jcycle = fp_isog_cycle((j,p),lk)
            qfcycle = qf_isog_cycle(qf,lk)
            assert len(jcycle) == len(qfcycle)
            for j1,qf1 in zip(jcycle,qfcycle):
                if j1 in j_to_qf:
                    assert j_to_qf[j1] == qf1
                else:
                    newassignments[j1] = qf1
        for j1,qf1 in newassignments.items():
            j_to_qf[j1] = qf1
            qf_to_j[qf1] = j1
    return j_to_qf


# Will replace previous with something built out of next function #

def ecqf_ord_n2_grdics(ap,ldata = None,ancdata = None):
    a,p = ap
    if a == 0:
        raise ValueError('Use ss alg')
    d = a*a-4*p
    if ldata == None:
        ldata = disc_rig_ssl_search(d)
        if type(ldata)!=dict or ldata['needs_sum']:
            raise ValueError('Need rig basis, use diff alg')
    ls = ldata['ls']
    ns = ldata['ns']
    if len(ns)>1 and max(ns[1:])>2:
        ls,ns = sortpair(ls,ns)
        if max(ns[1:])>2:
            raise ValueError('Multiple generators of order greater than 2')
    if ancdata == None:
        ancdata = get_ancestor_data_ord(ap)
    j0 = ancdata['leaves'][0]
    qf0 = class_group_id(d)
    if len(ls)==0:
        return {():(j0,qf0)}
    jcyc0 = fp_isog_cycle((j0,p),ls[0])
    qfcyc0 = qf_isog_cycle(qf0,ls[0])
    t_to_pairs = {}
    assert len(jcyc0)==len(qfcyc0)
    for i in range(len(jcyc0)):
        t_to_pairs[(i,)]=(jcyc0[i],qfcyc0[i])
    for l in ls[1:]:
        newdic = {}
        for t in t_to_pairs:
            jt,qft = t_to_pairs[t]
            jcyct =  fp_isog_cycle((jt,p),l)
            qfcyct = qf_isog_cycle(qft,l)
            for x in [0,1]:
                newdic[ext_tup(t,x)]=(jcyct[x],qfcyct[x])
        t_to_pairs =newdic
    return newdic

#####################
# (n,m) with n,m >2 #
# ###################

### Quadratic form (n,m) algorithms ###


def qf_oriented_cycle(qfs12:tuple[tuple[int]],l:int):
    qf1,qf2 = qfs12
    cycle = qf_isog_cycle(qf1,l)
    if len(cycle)<2:
        return []
    if cycle[1]==qf2:
        return cycle
    elif cycle[-1]==qf2:
        cycle_tail = cycle[1:]
        cycle = [cycle[0]]+cycle_tail[::-1]
        return cycle
    else:
        return []

def qf_l123_trios(qf0:tuple[int,int,int],l123:tuple[int,int,int])->list[tuple[int,int,int]]:
    l1,l2,l3 = l123
    qfs_l1 = qf_isogs_hor(qf0,l1)
    qfs_l2 = qf_isogs_hor(qf0,l2)
    qfs_l3 = qf_isogs_hor(qf0,l3)
    qftrios = []
    for qf1 in qfs_l1:
        qf1_l2_qfs = qf_isogs_hor(qf1,l2)
        for qf2 in qfs_l2:
            qfs_l123 = [qf for qf in qf_isogs_hor(qf2,l1) if qf in qf1_l2_qfs and qf in qfs_l3]
            if len(qfs_l123)>0:
                qftrios.append([qf1,qf2,qfs_l123[0]])
    return qftrios


def qf_l12_page(qf0,l12,qf123):
    l1,l2 = l12
    qf1,qf2,qf3 = tuple(qf123)
    row1 = qf_oriented_cycle((qf0,qf1),l1)
    row2 = qf_oriented_cycle((qf2,qf3),l1)
    n1 = len(row1)
    assert len(row2)==n1
    page = {}
    for i1 in range(n1):
        qfi1 = row1[i1]
        qfi2 = row2[i1]
        coli = qf_oriented_cycle((qfi1,qfi2),l2)
        for i2,qf in enumerate(coli):
            page[(i1,i2)] = qf
    return page

def qf_nn_page_from_ldata(d,ldata):
    if len(ldata['ls'])!= 2:
        raise ValueError('Use diff alg')
    l1, l2 = ldata['ls']
    if ldata['l_sum'] == None:
        raise ValueError('Use diff alg')
    l3 = ldata['l_sum']
    qf0 = class_group_id(d)
    trios_all = qf_l123_trios(qf0,(l1,l2,l3))
    if len(trios_all)==0:
        raise ValueError('No trios found')
    qf123 = trios_all[0]
    return qf_l12_page(qf0,(l1,l2),qf123)

### Fp (n,m) algorithms ###

def fp_oriented_cycle(j12,p,l):
    j1,j2 = j12
    cyc = fp_isog_cycle((j1,p),l)
    if len(cyc)<2:
        raise ValueError('No cycle')
    if cyc[1]==j2:
        return cyc
    elif cyc[-1]==j2:
        tail = cyc[1:]
        cyc = [cyc[0]]+tail[::-1]
        return cyc
    else:
        raise ValueError('Incompatible orientation')

def ecfp_l123_bases(j0,p,l123):
    l1,l2,l3 = l123
    js_l1 = fp_isog_codomains(j=j0,p=p,l=l1)
    js_l2 = fp_isog_codomains(j=j0,p=p,l=l2)
    js_l3 = fp_isog_codomains(j=j0,p=p,l=l3)
    jtrios = []
    for j1 in js_l1:
        j1_l2js = fp_isog_codomains(j=j1,p=p,l=l2)
        for j2 in js_l2:
            js_123 = [j for j in fp_isog_codomains(j=j2,p=p,l=l1) if j in j1_l2js and j in js_l3]
            if len(js_123)>0:
                jtrios.append([j1,j2,js_123[0]])
    return jtrios

def ecfp_l123_page(j0,p,l123):
    l1,l2,l3 = l123
    trios = ecfp_l123_bases(j0,p,l123)
    if len(trios)== 0:
        raise ValueError('No basis found')
    trio = trios[0]
    j1,j2,j3 = tuple(trio)
    row1 = fp_oriented_cycle((j0,j1),p,l1)
    row2 = fp_oriented_cycle((j2,j3),p,l1)
    n1 = len(row1)
    assert len(row2)==n1
    page = {}
    for i1 in range(n1):
        ji1 = row1[i1]
        ji2 = row2[i1]
        coli = fp_oriented_cycle((ji1,ji2),p,l2)
        for i2,j in enumerate(coli):
            page[(i1,i2)] = j
    return page

### Main Computation ###

def ecqf_ord_nn_grdics(ap,ldata = None,ancdata = None):
    a,p = ap
    d = a*a-4*p
    if a == 0:
        raise ValueError('Use supersingular alg')
    if ldata == None:
        ldata = disc_rig_ssl_search(d)
    if ldata['l_sum'] == None:
        raise ValueError('Need rig data for this alg, none available')
    ls = ldata['ls']
    if len(ls) != 2:
        raise ValueError('Use diff algo')
    l1,l2 = ls
    l3 = ldata['l_sum']
    if ancdata == None:
        ancdata = get_ancestor_data_ord(ap)
    j0 = ancdata['leaves'][0]
    j_dic = ecfp_l123_page(j0,p,(l1,l2,l3))
    qf_dic = qf_nn_page_from_ldata(d,ldata)
    return j_dic,qf_dic

def ecqf_ord_nn_bij_leaves(ap,ldata = None,ancdata = None):
    j_dic,qf_dic = ecqf_ord_nn_grdics(ap,ldata,ancdata)
    return compdiv_dics(j_dic,qf_dic)


########################
# End-to-end algorithm #
########################

# Given a bijection of leaves and ancestor data over Fp,
# this extends the bijection to the whole isogeny class

def vert_isog_ext(j_to_qf:dict,vertical_iso_data:dict)->dict:
    for l in vertical_iso_data['ancestor_data']:
        ancsl = vertical_iso_data['ancestor_data'][l]
        nextbatch = [j for j in j_to_qf if j in ancsl and ancsl[j] not in j_to_qf]
        while len(nextbatch)>0:
            currentbatch = nextbatch.copy()
            nextbatch = []
            for j0 in currentbatch:
                j1 = ancsl[j0]
                if j1 not in j_to_qf:
                    qf0 = j_to_qf[j0]
                    qf1 = qf_parents(qf0,l)[0]
                    j_to_qf[j1] = qf1
                    if j1 in ancsl:
                        nextbatch.append(j1)
    return j_to_qf



def ecqf_ord_bij(ap,ldata=None):
    a,p = ap
    if a == 0:
        raise ValueError('Use supersingular algorithm')
    d = a*a-4*p
    hd = small_bij_check(d)
    if len(hd)>0:
        return {j%p:hd[j] for j in hd}
    if ldata == None:
        ldata = disc_rig_ssl_search(d)
    if type(ldata)!=dict:
        raise ValueError('Failed to find ls with current tools')
    ancdata = get_ancestor_data_ord(ap)
    if not ldata['needs_sum']:
        leaf_bij = ecqf_ord_n2_bij_leaves(ap,vertical_iso_data=ancdata,ldata = ldata)
    elif len(ldata['ls'])==2 and ldata['l_sum']!= None:
        leaf_bij = ecqf_ord_nn_bij_leaves(ap,ldata=ldata,ancdata=ancdata)
    else:
        raise ValueError('Failed to compute leaf bijection with current tools')
    return vert_isog_ext(leaf_bij,ancdata)
    
