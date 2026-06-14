import itertools
import json
from pathlib import Path
from qfs import *
from nt import discfac,quad_rec,divisors,primesBetween
from identities import *
from ecfp import fp_isog_codomains, trfr_to_js
from modularpolynomials import *

_DATA_DIR = Path(__file__).parent / 'data'

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

with open(_DATA_DIR / 'ecqf_ord_pcbij_4to256.json', 'r') as f:
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


# ---------------------------------------------------------------------------
# The l-set search and the j <-> qf bijection computations that used to live
# below this point now have a single home in ecqf_bij.py:
#   * l-set search:          disc_rigid_lset_search
#   * end-to-end (a,p,ls):   ecqf_full_bijection_ord
#   * cache-aware (a,p):     rigid_cache.ecqf_ord_bij_cached(ap, use_cache=True)
# They were removed from here to drop the duplicate, divergent copies.
# ---------------------------------------------------------------------------
