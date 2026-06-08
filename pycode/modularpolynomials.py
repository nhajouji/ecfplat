## Data loading

import json
from pathlib import Path
from ringclasses import Polynomial

_DATA_DIR = Path(__file__).parent / 'data'

with open(_DATA_DIR / 'atkinpolys.json', 'r') as f:
    atkin_polys_raw = json.load(f)
atkin_polys_dict = {int(p):atkin_polys_raw[p] for p in atkin_polys_raw}

with open(_DATA_DIR / 'hilbpolys.json', 'r') as f:
    hilbpolys_raw = json.load(f)
hilb_polys_dict = {int(d):hilbpolys_raw[d] for d in hilbpolys_raw}
heeg_js = {d:-hilb_polys_dict[d][0] for d in hilb_polys_dict if len(hilb_polys_dict[d])==2}


## Polynomials and Rational Functions

def poly_eval_mod(coefs_lc_to_const:list[int],x:int,p:int,rev = False):
    if rev == True:
        coefs_lc_to_const = coefs_lc_to_const[::-1]
    evx = coefs_lc_to_const[0]%p
    for c in coefs_lc_to_const[1:]:
        evx = (evx*x+c)%p
    return evx


def rat_eval_mod(coefdic_lc_to_const:dict,x:int,p:int,rev = False):
    ev = 1
    for e in coefdic_lc_to_const:
        fac_coefs = coefdic_lc_to_const[e]
        fac_ev = poly_eval_mod(coefs_lc_to_const=fac_coefs,x=x,p=p,rev=rev)
        if e>= 0 or fac_ev %p != 0:
            ev = (ev*pow(fac_ev,e,p)) %p
        else:
            return 'Infinity'
    return ev


## Count roots

def count_roots_fp_bf(coefs_lc_to_const:list[int],p:int)->int:
    return len([x for x in range(p) if poly_eval_mod(coefs_lc_to_const,x,p)==0])

## Evaluating modular polynomials

def eval_atk(x:int,l:int,p:int):
    if l not in atkin_polys_dict:
        raise ValueError('l does not divide order of the Monster group')
    cs = atkin_polys_dict[l]
    a,b1,b3= tuple([poly_eval_mod(ci[::-1],x,p) for ci in cs])
    b = (b1*pow(b3,3,p))%p
    return ((-a)%p,b)


def atk_poly_a(l:int,p=0):
    if l not in atkin_polys_dict:
        raise ValueError('l does not divide order of the Monster group')
    return Polynomial(atkin_polys_dict[l][0],char=p)

#This, in isolation, captures endorphisms of degree l for j = 0
def atk_poly_b1(l:int,p=0):
    if l not in atkin_polys_dict:
        raise ValueError('l does not divide order of the Monster group')
    return Polynomial(atkin_polys_dict[l][1],char=p)

def atk_poly_b3(l:int,p=0):
    if l not in atkin_polys_dict:
        raise ValueError('l does not divide order of the Monster group')
    return Polynomial(atkin_polys_dict[l][2],char=p)

def atk_poly_b(l:int,p=0):
    b1 = atk_poly_b1(l,p)
    b3 = atk_poly_b3(l,p)
    return b1*(b3**3)

def atk_poly_ab(l:int,p=0):
    return (atk_poly_a(l,p),atk_poly_b(l,p))

def atk_at_j(j:int,l:int,p=0):
    j2 = Polynomial([j**2],char=p)
    a,b = atk_poly_ab(l,p)
    return j2-j*a+b


def atk_at_j_fpfac(j:int,l:int,p:int):
    poly = atk_at_j(j,l,p).mod(p)
    return poly.fp_factor()