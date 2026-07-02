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
from nt import primeQ

_DATA_DIR = Path(__file__).parent / 'data'
_EXT_PATH = _DATA_DIR / 'hilbert_roots_ext.json'


####################################
# Step 1a: harvest roots from data #
####################################

def class_roots_by_disc(j_to_qf: dict) -> dict:
    """{disc: sorted j's} for one class's j -> qf bijection."""
    groups = {}
    for j, qf in j_to_qf.items():
        groups.setdefault(qf_disc(qf), []).append(j)
    return {d: sorted(js) for d, js in groups.items()}


def cm_js_by_disc(source=None, include_ext: bool = True) -> dict:
    """{d: {p: sorted j's with endomorphism disc d over F_p}} from the bijection dicts.

    include_ext merges the beyond-1024 extension cache (hilbert_roots_ext.json).
    When two sources see the same (d, p) (twist classes, or base + extension), the
    root sets must agree; any mismatch means corrupted data, so it raises."""
    if source is None:
        source = ecqf_ord_1K_pc
    out = {}
    for (a, p), j_to_qf in source.items():
        if a <= 0:
            continue
        for d, js in class_roots_by_disc(j_to_qf).items():
            prev = out.setdefault(d, {}).setdefault(p, js)
            if prev != js:
                raise ValueError(f'Inconsistent roots for d={d}, p={p} (classes with traces +-{a})')
    if include_ext:
        for d, by_p in load_ext().items():
            for p, js in by_p.items():
                prev = out.setdefault(d, {}).setdefault(p, js)
                if prev != js:
                    raise ValueError(f'Extension cache disagrees with base data at d={d}, p={p}')
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


###############################################
# Extension beyond the p < 1024 precomputes   #
###############################################
# The precomputed bijections stop at p < 1024, which certifies only ~120 discs.
# For a target disc d we mint new CRT primes on demand: any prime p with
# a^2 - 4p = d*m^2 gives a class whose disc-d group is the full root set of
# H_d mod p.  The per-disc machinery (rigid l-set + qf labelling) is already
# cached in rigid_lset_cache.json down to -4096, so each new prime costs only
# the j-invariant side (ecqf_ord_bij_cached).  Results are persisted to
# hilbert_roots_ext.json as {d: {p: roots}}, including the groups of every
# other disc in the class closure -- those primes are free.

def load_ext(path=_EXT_PATH) -> dict:
    if not Path(path).exists():
        return {}
    with open(path) as f:
        raw = json.load(f)
    return {int(d): {int(p): js for p, js in by_p.items()} for d, by_p in raw.items()}


def save_ext(ext: dict, path=_EXT_PATH):
    tmp = str(path) + '.tmp'
    with open(tmp, 'w') as f:
        json.dump({str(d): {str(p): js for p, js in sorted(by_p.items())}
                   for d, by_p in sorted(ext.items(), reverse=True)}, f)
    Path(tmp).replace(path)


def crt_prime_candidates(d: int, pmin: int = 1024, pmax: int = 8192,
                         disc_floor: int = -4096) -> list[tuple]:
    """New CRT primes for d: [(p, a, D)] with a^2 - 4p = D = d*m^2, sorted by p.

    disc_floor keeps every class disc D inside the rigid-lset cache, so no new
    per-disc search is triggered.  Note m=1 can be structurally empty (d = 1 mod 8
    forces p even), which is why the m-loop matters."""
    out = []
    m = 1
    while d * m * m >= disc_floor:
        D = d * m * m
        a = int(math.isqrt(4 * pmin + D)) + 1 if 4 * pmin + D > 0 else 1
        if (a * a - D) % 2:
            a += 1                          # need a^2 = D (mod 2), then mod 4 is automatic
        while (p := (a * a - D) // 4) <= pmax:
            if p > pmin and primeQ(p):
                out.append((p, a, D))
            a += 2
        m += 1
    return sorted(out)


def harvest_class(a: int, p: int, ext: dict = None) -> dict:
    """Compute the (a, p) bijection and fold every disc group into ext."""
    from rigid_cache import ecqf_ord_bij_cached
    if ext is None:
        ext = {}
    groups = class_roots_by_disc(ecqf_ord_bij_cached((a, p)))
    for d, js in groups.items():
        prev = ext.setdefault(d, {}).setdefault(p, js)
        if prev != js:
            raise ValueError(f'Extension harvest inconsistent at d={d}, p={p}')
    return groups


def extend_disc(d: int, data: dict = None, pmax: int = 8192, save: bool = True,
                verbose: bool = True) -> dict:
    """Mint new primes for d until it is certified (or candidates run out).

    Returns the hilbert_via_crt record computed from the merged data.  Skips
    primes already known; records nothing when the class computation fails (rare
    rigid-l-set gaps) and moves on to the next candidate."""
    ext = load_ext()
    if data is None:
        data = cm_js_by_disc()              # already merges the current ext
    have = data.setdefault(d, {})
    bound = hilbert_bound_log2(d)
    bits = sum(math.log2(p) for p in have)
    for p, a, D in crt_prime_candidates(d, pmax=pmax):
        if bits > bound + 1:
            break
        if p in have:
            continue
        try:
            groups = harvest_class(a, p, ext)
        except (ValueError, KeyError) as e:
            if verbose:
                print(f'  skip (a={a}, p={p}, D={D}): {e}')
            continue
        for dd, js in groups.items():
            data.setdefault(dd, {}).setdefault(p, js)
        bits += math.log2(p)
        if verbose:
            print(f'  + p={p} (a={a}, D={D}): {bits:.1f} / {bound + 1:.1f} bits')
        if save:
            save_ext(ext)
    return hilbert_via_crt(d, data)
