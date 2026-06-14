from nt import primesBetween,discfac,quad_rec
from qfs import *
from ecfp import trfr_to_js
from graph_tools import compute_bijection_zn
from modularpolynomials import eval_atk,small_bij_check
from misctools import *
from ecqf import qf_reps_pm

import itertools
ssprimes = primesBetween(2,36)+[41]+[47+12*k for k in range(3)]



                ######################
                # Rigid l-set search #
                ######################


#################
# Search l-gens #
#################
ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]

###########################
# Obtaining generating ls #
###########################
ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]

def scoredic_to_tups(scores:dict):
    ls = [l for l in scores]
    ls.sort(reverse=True)
    return tuple(ls),tuple([scores[l] for l in ls])




def qf_reps_pm(d:int):
    return [qf for qf in get_qfs_strict(d) if qf[1]>=0]

def disc_to_ssls(d:int,ls=ssprimes):
    qf_reps = qf_reps_pm(d)
    ls = {}
    qf0 = class_group_id(d)
    for qf in qf_reps:
        qfls = [l for l in qf_evs_inrange(qf,10) if l in ssprimes]
        if len(qfls)>0:
            l0 = min(qfls)
            ls[l0]=len(qf_isog_cycle(qf0,l0))
    return ls

def disc_to_ssl_qfs(d,ls=ssprimes):
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

def score_tuple(ls,scores):
    prod = 1
    for l in ls:
        prod*=scores[l]
    return prod

def check_lset(d,ls):
    cld = len(get_qfs_strict(d))
    lgr_size = len(gen_qfs_d_ls(d,ls))
    prd = 1
    nlist = []
    for l in ls:
        nl = len(gen_qfs_d_ls(d,[l]))
        prd*=nl
    return [cld,lgr_size,prd]

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
    cld = len(get_qfs_strict(d))
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
    cld = len(get_qfs_strict(d))
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


                ####################################
                # Rigid l-set search (from scratch) #
                ####################################
#
# Goal: given the class group A of discriminant d and a pool of candidate
# primes `ls` (each prime l represents an element x_l in A via the horizontal
# isogeny action), find a subset whose elements form a *rigid spanning set*,
# packaged so it can be handed straight to `ecqf_full_bijection_ord` /
# `compute_bijection_zn`.
#
# What the downstream cube construction in graph_tools needs from the tuple it
# is given:
#   * The elements must be INDEPENDENT generators, i.e. a direct-product basis
#     A = prod <x_l>.  For a generating set this is equivalent to
#         prod_l ord(x_l) == |A|
#     (a generating set always over-counts; equality <=> direct product).
#     This is strictly stronger than "irredundant" but is what the labelling
#     algorithm requires, and it implies irredundancy.
#   * Generators of order 2 are handled separately (extend_bijection_zn2) and
#     just need to be independent 2-torsion; no orientation data is needed.
#   * If there are m >= 2 generators of order > 2, the cube needs one extra
#     "pinning" element to fix the relative orientation of the m cycles.
#     compute_bijection_zn uses the *last* long prime for this, and the tree
#     search only closes up if that prime's element is a signed sum
#         +- x_1 +- ... +- x_m
#     of the other long generators -- i.e. exactly the "sum of the order > 2
#     generators" x* from the notebook.  With m <= 1 no pinning element is
#     needed.
#
# So a valid rigid tuple is
#     ls_rig = (b_1, ..., b_m, l_sum, c_1, ..., c_k)
# with b_i the independent order > 2 generators, l_sum a prime representing a
# signed sum of them (omitted when m <= 1), and c_j independent order-2
# generators (their position is irrelevant -- l2_split pulls them out by order).


def _l_element_order(d, l, qf0=None):
    """Order of the element x_l = length of its horizontal isogeny cycle."""
    if qf0 is None:
        qf0 = class_group_id(d)
    return len(qf_isog_cycle(qf0, l))


def _rig_candidates(d, ls, qf0):
    """Map each usable prime to ord(x_l), dropping primes with trivial x_l."""
    cand = {}
    for l in ls:
        nontriv = [qf for qf in qf_isogs_hor(qf0, l) if qf != qf0]
        if len(nontriv) == 0:
            continue                       # x_l is the identity -> useless
        o = _l_element_order(d, l, qf0)
        if o > 1:
            cand[l] = o
    return cand


def _find_independent_basis(d, cand_sorted, cand, target):
    """DFS for an independent generating set (prod of orders == |A|).

    Adds a prime only when it extends the current subgroup *independently*:
        |<S + [l]>| == |<S>| * ord(x_l).
    Returns the chosen list of primes, or None if no such set exists in the
    pool.  cand_sorted should be ordered by descending order so the first
    basis found uses few, large-order generators (minimal rank)."""
    def dfs(start, S, g):
        if g == target:
            return list(S)
        for i in range(start, len(cand_sorted)):
            l = cand_sorted[i]
            g2 = len(gen_qfs_d_ls(d, S + [l]))
            if g2 == g * cand[l]:
                res = dfs(i + 1, S + [l], g2)
                if res is not None:
                    return res
        return None
    return dfs(0, [], 1)


def _best_independent_partial(d, cand_sorted, cand):
    """Largest independent set we can build -- used to report a near miss."""
    best = []
    best_g = 1
    def dfs(start, S, g):
        nonlocal best, best_g
        if g > best_g:
            best_g, best = g, list(S)
        for i in range(start, len(cand_sorted)):
            l = cand_sorted[i]
            g2 = len(gen_qfs_d_ls(d, S + [l]))
            if g2 == g * cand[l]:
                dfs(i + 1, S + [l], g2)
    dfs(0, [], 1)
    return best


def _signed_sum_corners(d, longs, qf0):
    """All cube far-corners  +- x_1 +- ... +- x_m  for the long generators.

    These are exactly the endpoints the cube tree-search produces, and the
    pinning prime l_sum must represent one of them."""
    frontier = {qf0}
    for l in longs:
        nxt = set()
        for qf in frontier:
            for qf1 in qf_isogs_hor(qf, l):     # qf +- x_l
                nxt.add(qf1)
        frontier = nxt
    return frontier


def disc_rigid_lset_search(d, ls=ssprimes):
    """Search for a rigid spanning l-set for the class group of discriminant d.

    Returns a dict.  On success (`success=True`) the key `ls_rig` is a flat
    tuple of primes ready to pass to `ecqf_full_bijection_ord(a, p, ls_rig)`:
        ls_rig = long basis ... , l_sum (if needed), order-2 generators ...
    Other keys: `ls_basis` (independent order>2 gens), `ls_2tors` (independent
    2-torsion gens), `l_sum`, `needs_sum`, `ns` (orders aligned with the
    generators, excluding the pinning prime), and `order` = |A|.

    On failure `success=False`, `message` says which step failed and `best`
    holds the best independent subset found."""
    qf0 = class_group_id(d)
    A = get_qfs_strict(d)
    order = len(A)
    out = {'d': d, 'order': order, 'success': False, 'message': '',
           'ls_basis': (), 'ls_2tors': (), 'l_sum': None, 'needs_sum': False,
           'ls': (), 'ns': (), 'ls_rig': ()}

    if order == 1:
        out.update(success=True, message='Trivial class group.')
        return out

    cand = _rig_candidates(d, ls, qf0)
    if len(cand) == 0:
        out['message'] = 'Failed (candidates): no prime in the pool represents a nontrivial element.'
        return out

    cand_sorted = sorted(cand, key=lambda l: cand[l], reverse=True)

    # Step 1: an independent generating set (a direct-product basis of A).
    basis = _find_independent_basis(d, cand_sorted, cand, order)
    if basis is None:
        best = _best_independent_partial(d, cand_sorted, cand)
        out.update(message='Failed (spanning): the pool has no independent '
                           'generating set for A.',
                   best=tuple(best),
                   best_order=len(gen_qfs_d_ls(d, list(best))))
        return out

    longs = sorted((l for l in basis if cand[l] > 2),
                   key=lambda l: cand[l], reverse=True)
    twos = [l for l in basis if cand[l] == 2]

    # Step 2: if >= 2 long generators, find a pinning prime for the sum element.
    l_sum = None
    if len(longs) >= 2:
        corners = _signed_sum_corners(d, longs, qf0)
        used = set(basis)
        for l in cand_sorted:
            if l in used:
                continue
            if any(qf in corners for qf in qf_isogs_hor(qf0, l) if qf != qf0):
                l_sum = l
                break
        if l_sum is None:
            out.update(message='Failed (rigid): found a basis with >= 2 '
                               'generators of order > 2, but no pool prime '
                               'represents the sum element x* = x_1 + ... + x_m.',
                       best=tuple(basis), ls_basis=tuple(longs),
                       ls_2tors=tuple(twos), needs_sum=True)
            return out

    # Step 3: assemble.  l_sum goes last among the long primes; 2-torsion last.
    ls_rig = tuple(longs) + ((l_sum,) if l_sum is not None else ()) + tuple(twos)
    out.update(success=True,
               message='Found a rigid spanning set.',
               ls_basis=tuple(longs), ls_2tors=tuple(twos),
               l_sum=l_sum, needs_sum=(l_sum is not None),
               ls=tuple(longs) + tuple(twos),
               ns=tuple(cand[l] for l in longs) + tuple(cand[l] for l in twos),
               ls_rig=ls_rig)
    return out


        #########################################
        # Step 2: Collect Relevant Isogeny Data #
        #########################################

# Quadratic form/lattice data

def qf_isog_data(d,ls):
    qfs = get_qfs_strict(d)
    # De-duplicate each neighbour list: for an order-2 element x_l the +x_l and
    # -x_l isogenies land on the same class, and qf_isogs_hor returns it twice.
    # Collapsing the duplicate lets l2_split correctly read such an l as a
    # length-1 (order-2) direction rather than mistaking it for a long cycle.
    return {l:{qf:list(dict.fromkeys(qf_isogs_hor(qf,l))) for qf in qfs} for l in ls}


######################
# ECFP Ordinary Data #
######################


def js_to_rabs(js,p):
    ab_to_js = {}
    for i,j1 in enumerate(js):
        for j2 in js[i:]:
            ab_to_js[((-j1-j2)%p,(j1*j2)%p)] = (j1,j2)
    return ab_to_js

def ecfp_nbr_data_ord_X1(ap,l,jdata = {}):
    a,p = ap
    if len(jdata) == 0:
        js = trfr_to_js(a,p)
        rabs = js_to_rabs(js,p)
    else:
        js = jdata['js']
        rabs = jdata['rabs']
    nbrdata = {j:[] for j in js}
    for x in range(p):
        evx = eval_atk(x,l,p)
        if evx in rabs:
            j1,j2 = rabs[evx]
            nbrdata[j1].append(j2)
            nbrdata[j2].append(j1)
    return nbrdata


def ecfp_nbr_data_ord(ap,ls):
    a,p = ap
    js = trfr_to_js(a,p)
    rabs = js_to_rabs(js,p)
    jspc = {'js':js,'rabs':rabs}
    # De-duplicate neighbour lists (an order-2 direction lists the same j twice),
    # mirroring qf_isog_data, so l2_split reads order-2 vs long the same way on the
    # j side as on the qf side.  The vertical/ancestor code keeps multiplicities by
    # calling ecfp_nbr_data_ord_X1 directly, so it is unaffected by this.
    out = {}
    for l in ls:
        nbr = ecfp_nbr_data_ord_X1(ap,l,jspc)
        out[l] = {j:list(dict.fromkeys(nbr[j])) for j in nbr}
    return out

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
        nbrs = ecfp_nbr_data_ord_X1(ap,2)
        for j in leaves:
            ancs[j]=nbrs[j][0]
        return {'ancestor_data':{2:ancs},'leaves':leaves,'js_all':js}
    ls = [l for l in primesBetween(1,c+1) if c%l ==0]
    anc_data = {}
    leaf_cands = [j for j in js if j!= 0 and (j-1728)%p !=0]
    for l in ls:
        nbrs_l = ecfp_nbr_data_ord_X1((a,p),l)
        anc_data[l]=tree_edges_to_ancestors(nbrs_l)
        leaf_cands = [j for j in leaf_cands if len(nbrs_l[j])==1]
    return {'ancestor_data':anc_data,'leaves':leaf_cands,'js_all':js}

def zn_ecqf_bij(a:int,p:int,ls:tuple[int]):
    anc_data = get_ancestor_data_ord((a,p))
    j0 = anc_data['leaves'][0]
    isog_data_horz_fp = ecfp_nbr_data_ord((a,p),ls)
    zn_to_j = compute_bijection_zn(ls,isog_data_horz_fp,j0)
    zn_to_qf = compute_bijection_zn(ls,qf_isog_data(a*a-4*p,ls),class_group_id(a*a-4*p))
    assert len(zn_to_j)==len(zn_to_qf)
    assert len({t for t in zn_to_qf if t not in zn_to_j})== 0
    return {t:(zn_to_j[t],zn_to_qf[t]) for t in zn_to_qf}


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

        #########################################
        # Orientation canonicalization          #
        #########################################
#
# A rigid l-set pins the cube labelling up to the single global automorphism
# x -> -x (complex conjugation / class-group inversion); see the rigid l-set
# notes above.  So each of the two labellings (qf side and j side) has exactly
# one binary orientation freedom.  We fix both with deterministic, intrinsic
# rules so the end-to-end bijection is reproducible and -- crucially -- the
# qf-side labelling is a function of d alone, hence consistent across every
# (a,p) in the same discriminant (and so across isogeny classes).
#
# Convention:
#   * qf side: the first long generator e_1 = (1,0,...,0) maps to the reduced
#     form with b > 0.  (e_1's image has order > 2, so it is non-ambiguous and
#     b != 0; its inverse is the b < 0 form.)
#   * j side: the root's +1 neighbour in coordinate 0 is the numerically
#     smaller j-invariant.
# When there is no generator of order > 2 the group is 2-torsion, every form is
# ambiguous, the j-invariants are all "real", and conjugation is trivial -- both
# rules then leave the labelling untouched.

def _labelling_orders(bij):
    keys = list(bij)
    n = len(keys[0])
    return tuple(1 + max(k[i] for k in keys) for i in range(n))

def _neg_tuple(t,orders):
    return tuple((-ti) % n for ti,n in zip(t,orders))

def canonicalize_qf_labelling(zn_to_qf):
    """Pin the global x->-x freedom on the qf side (b>0 at the first long gen)."""
    if not zn_to_qf:
        return zn_to_qf
    n = len(next(iter(zn_to_qf)))
    if n == 0:
        return zn_to_qf
    e1 = (1,)+(0,)*(n-1)
    q1 = zn_to_qf.get(e1)
    if q1 is None or class_group_inv(q1) == q1:   # ambiguous => no order>2 coord
        return zn_to_qf
    if q1[1] < 0:
        return {t:class_group_inv(v) for t,v in zn_to_qf.items()}
    return zn_to_qf

def canonicalize_j_labelling(zn_to_j):
    """Pin the same global freedom on the j side (smaller j at root's +1 nbr)."""
    if not zn_to_j:
        return zn_to_j
    orders = _labelling_orders(zn_to_j)
    if len(orders) == 0 or orders[0] <= 2:        # coord 0 is 2-torsion => trivial
        return zn_to_j
    n = len(orders)
    e1  = (1,)+(0,)*(n-1)
    em1 = (orders[0]-1,)+(0,)*(n-1)
    if zn_to_j[e1] > zn_to_j[em1]:
        return {_neg_tuple(t,orders):v for t,v in zn_to_j.items()}
    return zn_to_j


def ecqf_full_bijection_ord(a:int,p:int,ls:tuple[int],zn_to_qf=None):
    # zn_to_qf is the (a,p)-independent labelling tuple -> qf of the class group.
    # It depends only on (d, ls), so it can be precomputed once per discriminant
    # and reused for every (a,p); pass it in to skip the expensive recompute.
    if a == 0:
        raise ValueError('Use supersingular algorithm instead')
    d = a*a-4*p
    hd = small_bij_check(d)
    if len(hd)>0:
        return {j%p:hd[j] for j in hd}
    vert_iso_data_fp = get_ancestor_data_ord((a,p))
    j0 = vert_iso_data_fp['leaves'][0]
    isog_data_horz_fp = ecfp_nbr_data_ord((a,p),ls)
    zn_to_j = canonicalize_j_labelling(compute_bijection_zn(ls,isog_data_horz_fp,j0))
    if zn_to_qf is None:
        zn_to_qf = compute_bijection_zn(ls,qf_isog_data(d,ls),class_group_id(d))
    zn_to_qf = canonicalize_qf_labelling(zn_to_qf)   # idempotent if already canonical
    assert len(zn_to_j)==len(zn_to_qf)
    assert len({t for t in zn_to_qf if t not in zn_to_j})== 0
    return vert_isog_ext(j_to_qf=compdiv_dics(zn_to_j,zn_to_qf),vertical_iso_data=vert_iso_data_fp)

