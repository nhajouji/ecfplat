"""Hilbert class polynomials H_d(x) via CRT from the precomputed bijection data.

Steps 1-2 of the CRT modular-polynomial project.  For each ordinary class (a, p)
in ecqf_ord_1K_pc, the bijection assigns every j-invariant a form whose
discriminant is the discriminant of the curve's endomorphism order; the j's
carrying discriminant d are exactly the roots of H_d mod p.  So

    H_d(x)  =  prod (x - j)   (mod p)

for every prime p with a^2 - 4p = d * m^2 for some a, m, and CRT across those
primes recovers the integer coefficients once the modulus exceeds twice the
largest coefficient.  The "enough primes" criterion uses the numeric bound

    max |coef of H_d| <= prod over reduced forms (1 + |j(tau)|),

with j(tau) approximated from the Fourier expansion (data/jq_coeffs.json).  The
Fourier series enters only through this bound, never the exact computation.
"""

import math
import cmath
import json
from pathlib import Path

from alg_classes import GF_p, poly_ring, poly_crt, Poly
from ecqf_tools import ecqf_ord_1K_pc
from qfs import qf_disc, get_qfs_strict
from modularpolynomials import hilb_polys_dict, _jc

_DATA_DIR = Path(__file__).parent / 'data'


####################################
# Step 1a: harvest roots from data #
####################################

def cm_js_by_disc(source=None) -> dict:
    """{d: {p: sorted j's with endomorphism disc d over F_p}} from the bijection dicts.

    When two classes (a1, p), (a2, p) both see disc d (twists), the root sets must
    agree; any mismatch means corrupted data, so it raises."""
    if source is None:
        source = ecqf_ord_1K_pc
    out = {}
    for (a, p), j_to_qf in source.items():
        if a <= 0:
            continue
        groups = {}
        for j, qf in j_to_qf.items():
            groups.setdefault(qf_disc(qf), []).append(j)
        for d, js in groups.items():
            js = sorted(js)
            prev = out.setdefault(d, {}).setdefault(p, js)
            if prev != js:
                raise ValueError(f'Inconsistent roots for d={d}, p={p} (classes with traces +-{a})')
    return out


##################################
# Step 1b: H_d mod p from roots  #
##################################

def hilbert_mod_p(js: list[int], p: int) -> Poly:
    """prod (x - j) over F_p, as a Poly over GF_p(p)."""
    return poly_ring(GF_p(p)).from_roots(js)


######################################
# Step 1c: coefficient bound (float) #
######################################

def j_numeric(tau: complex, nmax: int = 400, tol_log: float = -40.0) -> complex:
    """j(tau) from the Fourier expansion; adequate for bounds, not exact work."""
    q = cmath.exp(2j * cmath.pi * tau)
    lq = math.log(abs(q))
    j = 1 / q + 744
    for n in range(1, nmax + 1):
        c = _jc(n)
        if math.log(c) + n * lq < tol_log:      # term magnitude in logs (c(n) > 0)
            break
        j += c * q ** n
    return j


def cm_points(d: int) -> list[complex]:
    """The h(d) CM points tau = (-b + sqrt(d))/2a, one per primitive reduced form."""
    sq = math.sqrt(abs(d))
    return [complex(-b, sq) / (2 * a) for a, b, c in get_qfs_strict(d)]


def hilbert_bound_log2(d: int, safety_bits: float = 8.0) -> float:
    """log2 upper bound on max |coefficient of H_d|.

    Coefficients are elementary symmetric functions of the CM j-values, so each is
    at most prod (1 + |j_i|); summing log2(1 + |j_i|) avoids float overflow.  The
    safety bits absorb the truncation/rounding error of the numeric j's."""
    return sum(math.log2(1 + abs(j_numeric(tau))) for tau in cm_points(d)) + safety_bits


#########################
# Step 1d: CRT assembly #
#########################

def hilbert_via_crt(d: int, data: dict = None, use_all_primes: bool = False) -> dict:
    """Assemble H_d by CRT from the harvested data.

    Returns a record with the balanced-CRT coefficients (low-to-high), the primes
    used, and certified=True when the modulus provably exceeds twice the
    coefficient bound -- i.e. the coefficients are the true integers."""
    if data is None:
        data = cm_js_by_disc()
    if d not in data:
        return {'d': d, 'status': 'no data'}
    by_p = data[d]
    h = len(get_qfs_strict(d))
    for p, js in by_p.items():
        if len(js) != h:
            raise ValueError(f'd={d}, p={p}: {len(js)} roots but h(d)={h}')
    bound = hilbert_bound_log2(d)
    primes, bits = [], 0.0
    for p in sorted(by_p, reverse=True):        # largest primes first
        primes.append(p)
        bits += math.log2(p)
        if bits > bound + 1 and not use_all_primes:
            break
    F, M = poly_crt([(hilbert_mod_p(by_p[p], p), p) for p in primes])
    return {'d': d, 'status': 'ok', 'h': h, 'coefs': F.int_coefs(),
            'primes': primes, 'primes_available': len(by_p),
            'log2_modulus': bits, 'log2_bound': bound,
            'certified': bits > bound + 1}


##########################################
# Step 2: build / validate the library   #
##########################################

def certification_report(data: dict = None) -> dict:
    """Cheap pass over all harvested discs: how many primes we have vs how many we
    need.  No CRT is performed; use this to see what the current data can certify."""
    if data is None:
        data = cm_js_by_disc()
    report = {}
    for d in sorted(data, reverse=True):
        bound = hilbert_bound_log2(d)
        bits = sum(math.log2(p) for p in data[d])
        report[d] = {'h': len(get_qfs_strict(d)), 'primes': len(data[d]),
                     'log2_available': bits, 'log2_bound': bound,
                     'certified': bits > bound + 1}
    return report


def check_against_known(data: dict = None, verbose: bool = True) -> dict:
    """Rebuild every disc in hilbpolys.json that the data certifies; compare."""
    if data is None:
        data = cm_js_by_disc()
    results = {'match': [], 'mismatch': [], 'uncertified': [], 'no_data': []}
    for d in sorted(hilb_polys_dict, reverse=True):
        rec = hilbert_via_crt(d, data)
        if rec['status'] != 'ok':
            results['no_data'].append(d)
        elif not rec['certified']:
            results['uncertified'].append(d)
        elif rec['coefs'] == hilb_polys_dict[d]:
            results['match'].append(d)
        else:
            results['mismatch'].append(d)
    if verbose:
        print({k: len(v) for k, v in results.items()})
        if results['mismatch']:
            print('MISMATCHES:', results['mismatch'])
    return results


def build_hilbert_library(data: dict = None, certified_only: bool = True) -> dict:
    """{d: coefficients of H_d, low-to-high} for every disc the data supports."""
    if data is None:
        data = cm_js_by_disc()
    lib = {}
    for d in sorted(data, reverse=True):
        rec = hilbert_via_crt(d, data)
        if rec['status'] == 'ok' and (rec['certified'] or not certified_only):
            lib[d] = rec['coefs']
    return lib


def save_hilbert_library(lib: dict, path=None):
    """Write the library alongside hilbpolys.json (same format: {str(d): coefs})."""
    if path is None:
        path = _DATA_DIR / 'hilbpolys_crt.json'
    with open(path, 'w') as f:
        json.dump({str(d): lib[d] for d in sorted(lib, reverse=True)}, f)
    return path
