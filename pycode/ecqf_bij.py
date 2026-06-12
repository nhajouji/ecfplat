from nt import primesBetween,discfac,quad_rec
from qfs import class_group_id,get_qfs_strict,qf_isogs_hor
from ecfp import trfr_to_js
from graph_tools import compute_bijection_zn
from modularpolynomials import eval_atk,small_bij_check
from qfs import qf_parents
from misctools import *

ssprimes = primesBetween(2,36)+[41]+[47+12*k for k in range(3)]



                ######################
                # Rigid l-set search #
                ######################

from nt import primesBetween, discfac, quad_rec
from qfs import get_qfs_all, qf_isogs,qf_isogs_hor,qf_parents,qf_isog_cycle,class_group_id,qf_evs_inrange
from ecqf import qf_reps_pm

def qfs_in_conn_comp(qf0,ls):
    qfs = [qf0]
    n_ls = 1
    for l in ls:
        n_l = len(qf_isog_cycle(qf0,l))
        if n_l>1:
            n_ls*= n_l
            qfs_new = []
            for qf in qfs:
                qfs_new+=qf_isog_cycle(qf,l)
            qfs = list(set(qfs_new))
    return qfs,n_ls

def qfd_lgr_sizes(d,ls):
    qf0 = class_group_id(d)
    qfs,n = qfs_in_conn_comp(qf0,ls)
    return len(qfs),n

def qf_mimimize_ls(d,ls):
    gr_s1, gr_s2 = qfd_lgr_sizes(d,ls)
    lls= list(ls)
    i = 0
    while i < len(lls) and gr_s1 != gr_s2:
        lls_0 = lls[:i]+lls[i+1:]
        gr_s10,gr_s20 = qfd_lgr_sizes(d,tuple(lls_0))
        if gr_s10 == gr_s1:
            lls = lls_0
            gr_s1 = gr_s10
            gr_s2 = gr_s20
        else:
            i+=1
    return ls

        
def disc_to_ssls_all(d):
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

def disc_to_ssls_gens(d):
    return qf_mimimize_ls(d,[l for l in disc_to_ssls_all(d)])


        #########################################
        # Step 2: Collect Relevant Isogeny Data #
        #########################################

# Quadratic form/lattice data

def qf_isog_data(d,ls):
    qfs = get_qfs_strict(d)
    return {l:{qf:qf_isogs_hor(qf,l) for qf in qfs} for l in ls}


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
    return {l:ecfp_nbr_data_ord_X1(ap,l,jspc) for l in ls}

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

def ecqf_full_bijection_ord(a:int,p:int,ls:tuple[int]):
    if a == 0:
        raise ValueError('Use supersingular algorithm instead')
    d = a*a-4*p
    hd = small_bij_check(d)
    if len(hd)>0:
        return {j%p:hd[j] for j in hd}
    vert_iso_data_fp = get_ancestor_data_ord((a,p))
    j0 = vert_iso_data_fp['leaves'][0]
    isog_data_horz_fp = ecfp_nbr_data_ord((a,p),ls)
    zn_to_j = compute_bijection_zn(ls,isog_data_horz_fp,j0)
    zn_to_qf = compute_bijection_zn(ls,qf_isog_data(a*a-4*p,ls),class_group_id(a*a-4*p))
    assert len(zn_to_j)==len(zn_to_qf)
    assert len({t for t in zn_to_qf if t not in zn_to_j})== 0
    return vert_isog_ext(j_to_qf=compdiv_dics(zn_to_j,zn_to_qf),vertical_iso_data=vert_iso_data_fp)

