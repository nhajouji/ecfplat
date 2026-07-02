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
from nt import primeQ, discfac, primefact
from ecqf_bij import get_ancestor_data_ord, ssprimes

ATKIN_SET = frozenset(ssprimes)     # the 15 primes with an Atkin modular polynomial

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
# H_d mod p.  Identification of that group does NOT need the class-group
# labelling (rigid l-sets): the ancestor data of ecqf_bij determines every
# curve's endomorphism disc from the volcano structure alone, needing
# l-isogeny graphs only for the primes l dividing the conductor.  Results are
# persisted to hilbert_roots_ext.json as {d: {p: roots}}, including the groups
# of every other disc in the class closure -- those primes are free.

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


                # Criterion 1: which primes see disc d? #

def find_aps(d: int, N: int, pmin: int = 1024) -> list[tuple]:
    """All (p, a, m) with p prime in (pmin, N], a > 0 and a^2 - 4p = m^2 * d --
    i.e. the primes where F_p has ordinary curves with endomorphism disc d
    (d lies in the discriminant closure of a^2 - 4p exactly when such an m
    exists).  Sorted by p.  Note m=1 can be structurally empty (d = 1 mod 8
    forces p even), which is why the m-loop matters."""
    if d >= 0 or d % 4 > 1:
        raise ValueError(f'{d} is not a negative discriminant')
    out = []
    m = 1
    while m * m * abs(d) < 4 * N:
        D = d * m * m
        a = math.isqrt(4 * pmin + D) + 1 if 4 * pmin + D > 0 else 1
        if (a * a - D) % 2:
            a += 1                          # need a^2 = D (mod 2); mod 4 is then automatic
        while (p := (a * a - D) // 4) <= N:
            if p > pmin and primeQ(p):
                out.append((p, a, m))
            a += 2
        m += 1
    return sorted(out)


def crt_prime_candidates(d: int, pmin: int = 1024, pmax: int = 8192,
                         disc_floor: int = -4096) -> list[tuple]:
    """Legacy enumeration [(p, a, D)] kept for compatibility; use find_aps."""
    return [(p, a, d * m * m) for p, a, m in find_aps(d, pmax, pmin)
            if d * m * m >= disc_floor]


                # Criterion 2: can we identify the curves? #

def ap_status(a: int, p: int) -> dict:
    """Can every curve's endomorphism ring in class (a, p) be identified?

    The ancestor data needs the in-class l-isogeny graph for each prime l
    dividing the conductor of a^2 - 4p.  An Atkin prime (the 15 with a modular
    polynomial) supports any volcano depth.  For other l the only candidate is
    Velu's ramified/ascending branch, which is not implemented yet -- flagged
    fixable when the depth is 1; depth >= 2 needs kernel points over extensions
    of degree ~ l*(l-1), which is out of reach regardless.

    TODO (subring elimination, no isogenies needed): when the target order is
    identified by discarding the OTHER endo discs of the class, known Hilbert
    polynomials substitute for the l-graph -- e.g. conductor 37 has only two
    possible endo discs, and the disc-d0 curves are the roots of H_d0 mod p, so
    the disc-(d0*37^2) curves are the complement.  Works whenever H_d' is known
    for every proper superorder disc d' of the target."""
    d0, c = discfac(a * a - 4 * p)
    for l, e in sorted(primefact(c).items()):
        if l in ATKIN_SET:
            continue
        if e == 1:
            return {'ok': False, 'fixable': True,
                    'why': f'needs the ramified Velu ascent for l={l} || conductor (not implemented)'}
        return {'ok': False, 'fixable': False,
                'why': f'{l}^2 | conductor and no Atkin polynomial for {l}'}
    return {'ok': True, 'fixable': True, 'why': ''}


def class_endo_discs(a: int, p: int) -> dict:
    """{endo disc: sorted j's} across the whole ordinary class (a, p).

    Identification via ancestor data alone: the leaves of the volcano carry
    disc a^2 - 4p, and each ancestor step divides the disc by l^2.  This needs
    no rigid l-set and no class-group labelling, so it also works on class
    discs whose full bijection is blocked.  Every group is validated against
    its class number before returning."""
    D = a * a - 4 * p
    anc = get_ancestor_data_ord((a, p))
    amap = anc['ancestor_data']
    if amap and not isinstance(next(iter(amap.values())), dict):
        amap = {2: amap}                    # the (d0, c) = (-4, 2) branch returns a bare
                                            # {child: parent} map, not {l: {child: parent}}
    disc_of = {j: D for j in anc['leaves']}
    for l, ancsl in amap.items():           # ascend one prime at a time (as vert_isog_ext)
        batch = [j for j in disc_of if j in ancsl and ancsl[j] not in disc_of]
        while batch:
            nxt = []
            for j0 in batch:
                j1 = ancsl[j0]
                if j1 not in disc_of:
                    disc_of[j1] = disc_of[j0] // (l * l)
                    if j1 in ancsl:
                        nxt.append(j1)
            batch = nxt
    missing = set(anc['js_all']) - set(disc_of)
    if missing:
        raise ValueError(f'class ({a},{p}): {len(missing)} curves not reached by ancestor data')
    groups = {}
    for j, dd in disc_of.items():
        groups.setdefault(dd, []).append(j)
    for dd, js in groups.items():
        if len(js) != len(get_qfs_strict(dd)):
            raise ValueError(f'class ({a},{p}): disc {dd} has {len(js)} curves, h(d) = '
                             f'{len(get_qfs_strict(dd))}')
    return {dd: sorted(js) for dd, js in groups.items()}


def harvest_class(a: int, p: int, ext: dict = None) -> dict:
    """Identify the endo discs across class (a, p) and fold every group into ext."""
    if ext is None:
        ext = {}
    groups = class_endo_discs(a, p)
    for d, js in groups.items():
        prev = ext.setdefault(d, {}).setdefault(p, js)
        if prev != js:
            raise ValueError(f'Extension harvest inconsistent at d={d}, p={p}')
    return groups


                # Search report and end-to-end computation #

def hilbert_search_report(d: int, N: int, pmin: int = 1024, data: dict = None) -> dict:
    """Both criteria plus the bits budget: what does certifying H_d via primes
    <= N take, and does the search supply it?  No class is computed here."""
    if data is None:
        data = cm_js_by_disc()
    have = data.get(d, {})
    classes = [{'p': p, 'a': a, 'm': m, **ap_status(a, p)} for p, a, m in find_aps(d, N, pmin)]
    new_ps = sorted({c['p'] for c in classes if c['ok'] and c['p'] not in have})
    bound = hilbert_bound_log2(d)
    bits_have = sum(math.log2(p) for p in have)
    bits_new = sum(math.log2(p) for p in new_ps)
    can = bits_have + bits_new > bound + 1
    blocked = [c for c in classes if not c['ok']]
    verdict = (f'ok: {bits_have:.0f} bits on hand + {bits_new:.0f} available from '
               f'{len(new_ps)} new primes vs {bound + 1:.0f} needed' if can else
               f'FAIL: {bits_have:.0f} bits on hand + {bits_new:.0f} available < '
               f'{bound + 1:.0f} needed; raise N' +
               (f' (or unblock {len(blocked)} classes: {blocked[0]["why"]})' if blocked else ''))
    return {'d': d, 'N': N, 'h': len(get_qfs_strict(d)), 'log2_bound': bound,
            'primes_have': sorted(have), 'bits_have': bits_have,
            'new_usable_primes': new_ps, 'bits_new': bits_new,
            'classes': classes, 'blocked': blocked, 'can_certify': can, 'verdict': verdict}


def hilbert_poly_search(d: int, N: int, pmin: int = 1024, data: dict = None,
                        save: bool = True, verbose: bool = True, partial: bool = False) -> dict:
    """End-to-end: decide feasibility for H_d from primes in (pmin, N], then find
    the disc-d curves prime by prime and CRT the Hilbert polynomial.

    Returns the (certified) hilbert_via_crt record, or {'status': 'fail', ...}
    carrying the search report when N is too small.  partial=True harvests
    whatever is usable even when certification is out of reach."""
    if data is None:
        data = cm_js_by_disc()
    rep = hilbert_search_report(d, N, pmin, data)
    if not rep['can_certify'] and not partial:
        return {'status': 'fail', 'reason': rep['verdict'], 'report': rep}
    ext = load_ext()
    have = data.setdefault(d, {})
    bits, target = rep['bits_have'], rep['log2_bound'] + 1
    for cl in rep['classes']:               # cheapest primes first
        if bits > target:
            break
        p, a = cl['p'], cl['a']
        if not cl['ok'] or p in have:
            continue
        try:
            groups = harvest_class(a, p, ext)
        except ValueError as e:
            if verbose:
                print(f'  skip (a={a}, p={p}): {e}')
            continue
        for dd, js in groups.items():
            data.setdefault(dd, {}).setdefault(p, js)
        bits += math.log2(p)
        if verbose:
            print(f'  + p={p} (a={a}, m={cl["m"]}): {bits:.1f} / {target:.1f} bits')
        if save:
            save_ext(ext)
    return hilbert_via_crt(d, data)


def extend_disc(d: int, data: dict = None, pmax: int = 8192, save: bool = True,
                verbose: bool = True) -> dict:
    """Legacy entry point: harvest whatever primes <= pmax allow, certified or not."""
    return hilbert_poly_search(d, pmax, data=data, save=save, verbose=verbose, partial=True)
