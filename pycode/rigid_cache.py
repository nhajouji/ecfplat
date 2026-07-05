"""
Precomputed cache of rigid l-set search data and the qf-side bijection.

Motivation
----------
For a fixed discriminant d the expensive, class-group-side computations
    * the rigid l-set search        (disc_rigid_lset_search)
    * the labelling tuple -> qf      (zn_to_qf = compute_bijection_zn over qfs)
depend ONLY on d, not on the particular pair (a, p) with a^2 - 4p = d.  Both are
fully deterministic, so they can be precomputed once per discriminant, stored in
JSON, and reused for every (a, p) of that discriminant.  Only the j-invariant
side of the bijection (which does depend on (a, p)) is then recomputed per call.

JSON layout (data/rigid_lset_cache.json)
----------------------------------------
{
  "meta": {"pool": [...], "include_bijection": true, "updated": "..."},
  "discriminants": {
      "-24300": {
          "order": 54, "success": true, "message": "...",
          "ls_rig": [7,19,31,13], "ls_basis": [7,19,31], "ls_2tors": [],
          "l_sum": 13, "needs_sum": true, "ns": [6,3,3],
          "zn_to_qf": [[[x1,...,xn],[a,b,c]], ...]   # optional
      },
      ...
  }
}

CLI
---
    python rigid_cache.py --min -4000               # populate -3 .. -4000
    python rigid_cache.py --min -4000 --no-bijection # search data only (smaller)
    python rigid_cache.py --min -8000 --force        # recompute everything
Re-running only fills in what is missing, so it doubles as the "update later"
script: widen --min (or pass --force) when you have more compute.
"""

import os
import json
import datetime

from qfs import class_group_id, get_qfs_strict, qf_isogs_hor, qf_isog_cycle
from graph_tools import compute_bijection_zn
from modularpolynomials import small_bij_check, modular_prime_pool
from ecqf_bij import (disc_rigid_lset_search, qf_isog_data, ssprimes,
                      ecqf_full_bijection_ord, canonicalize_qf_labelling)

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
DEFAULT_CACHE = os.path.join(DATA_DIR, 'rigid_lset_cache.json')


def _ls_descs(ls):
    """Restore generator descriptors from JSON: lists -> tuples, recursively
    (powers [l,k], free generators ['sib',q]/['lift',q], and relaxed pins
    ['pin', pdesc, pairs, taus] with nested structure)."""
    def conv(x):
        if isinstance(x, list):
            return tuple(conv(y) for y in x)
        return x
    return tuple(conv(x) for x in ls)


#################
# Serialization #
#################

def _bij_to_json(bij):
    """zn_to_qf {tuple: (a,b,c)} -> JSON-friendly list of [key, qf] pairs."""
    return [[list(t), list(qf)] for t, qf in bij.items()]


def _bij_from_json(pairs):
    """Inverse of _bij_to_json."""
    return {tuple(t): tuple(qf) for t, qf in pairs}


#######################
# Compute one d-entry #
#######################

def compute_disc_entry(d, pool=None, include_bijection=True):
    """Run the search (and optionally the qf-side bijection) for one d.
    pool defaults to modular_prime_pool() (Atkin + cached classical Phi_l)."""
    res = disc_rigid_lset_search(d, pool)
    entry = {
        'order': res['order'],
        'success': res['success'],
        'message': res['message'],
        'ls_rig': list(res.get('ls_rig', ())),
        'ls_basis': list(res.get('ls_basis', ())),
        'ls_2tors': list(res.get('ls_2tors', ())),
        'l_sum': res.get('l_sum'),
        'needs_sum': res.get('needs_sum', False),
        'ns': list(res.get('ns', ())),
    }
    if not res['success']:
        entry['best'] = list(res.get('best', ()))
        return entry
    if include_bijection:
        ls = res['ls_rig']
        bij = compute_bijection_zn(ls, qf_isog_data(d, ls), class_group_id(d))
        bij = canonicalize_qf_labelling(bij)        # d-intrinsic orientation
        entry['zn_to_qf'] = _bij_to_json(bij)
    return entry


#############
# Load/save #
#############

_LOADED = {}                            # path -> parsed cache; parsing the (large and
                                        # growing) JSON per lookup dominated whole batches

def load_cache(path=DEFAULT_CACHE):
    if path in _LOADED:
        return _LOADED[path]
    if os.path.exists(path):
        with open(path) as f:
            cache = json.load(f)
    else:
        cache = {}
    cache.setdefault('meta', {})
    cache.setdefault('discriminants', {})
    _LOADED[path] = cache
    return cache


def save_cache(cache, path=DEFAULT_CACHE):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(cache, f)
    os.replace(tmp, path)               # atomic on POSIX
    _LOADED[path] = cache


###########
# Lookups #
###########

def get_disc_entry(d, cache=None, path=DEFAULT_CACHE, compute_if_missing=True,
                   persist=False, pool=None, include_bijection=True):
    """Return the cached entry for d, computing (and optionally storing) it if absent."""
    own = cache is None
    if own:
        cache = load_cache(path)
    discs = cache['discriminants']
    key = str(d)
    if key in discs:
        return discs[key]
    if not compute_if_missing:
        return None
    entry = compute_disc_entry(d, pool, include_bijection)
    if persist:
        discs[key] = entry
        save_cache(cache, path)
    return entry


def get_zn_to_qf(d, cache=None, path=DEFAULT_CACHE, **kw):
    """The labelling tuple -> qf for d, from cache or (re)computed from ls_rig."""
    entry = get_disc_entry(d, cache, path=path, **kw)
    if entry is None or not entry.get('success'):
        return None
    if 'zn_to_qf' in entry:
        return _bij_from_json(entry['zn_to_qf'])
    ls = _ls_descs(entry['ls_rig'])
    return canonicalize_qf_labelling(
        compute_bijection_zn(ls, qf_isog_data(d, ls), class_group_id(d)))


########################
# End-to-end (cached)  #
########################

def ecqf_ord_bij_cached(ap, use_cache=True, cache=None, path=DEFAULT_CACHE,
                        compute_if_missing=True):
    """Ordinary j -> qf bijection for ap=(a,p) -- the replacement for the old
    ecqf.ecqf_ord_bij, built on ecqf_bij's pipeline.

    use_cache=True (default): reuse the per-discriminant cache, so only the
    (a, p)-specific j-invariant side is computed here; falls back to an on-the-fly
    search when d is absent (unless compute_if_missing=False).
    use_cache=False: ignore the JSON cache entirely and compute live every call
    (disc_rigid_lset_search + ecqf_full_bijection_ord).
    Raises a clear error if no rigid l-set exists for d in the current prime pool."""
    a, p = ap
    if a == 0:
        raise ValueError('Use the supersingular algorithm for a = 0')
    d = a * a - 4 * p
    hd = small_bij_check(d)
    if len(hd) > 0:                       # tiny class number: closed-form table
        return {j % p: hd[j] for j in hd}
    if not use_cache:
        res = disc_rigid_lset_search(d)
        if not res['success']:
            raise ValueError(f"No rigid l-set for d = {d}: {res['message']}")
        return ecqf_full_bijection_ord(a, p, res['ls_rig'])
    entry = get_disc_entry(d, cache, path=path, compute_if_missing=compute_if_missing)
    if entry is None:
        raise KeyError(f'd = {d} is not in the cache and compute_if_missing=False')
    if not entry['success']:
        raise ValueError(f"No rigid l-set for d = {d}: {entry['message']}")
    ls = _ls_descs(entry['ls_rig'])
    zn_to_qf = _bij_from_json(entry['zn_to_qf']) if 'zn_to_qf' in entry else None
    return ecqf_full_bijection_ord(a, p, ls, zn_to_qf=zn_to_qf)


#######################
# Populate / update   #
#######################

def discriminants_in_range(dmin, dmax=-3):
    """Negative discriminants (== 0,1 mod 4) from dmax down to dmin, inclusive."""
    return [d for d in range(dmax, dmin - 1, -1) if d % 4 in (0, 1)]


def _entry_is_complete(entry, include_bijection):
    if not entry.get('success'):
        return True                       # a recorded failure is "done"
    if include_bijection:
        return 'zn_to_qf' in entry
    return True


def populate(dmin, dmax=-3, path=DEFAULT_CACHE, pool=None,
             include_bijection=True, force=False, save_every=200, verbose=True):
    """Fill the cache for every discriminant in [dmin, dmax].

    Skips discriminants already complete (unless force=True), so re-running with a
    wider range or after a speed-up only computes what is missing.  pool defaults
    to modular_prime_pool() at call time (Atkin + cached classical Phi_l)."""
    if pool is None:
        pool = modular_prime_pool()
    cache = load_cache(path)
    cache['meta'].update(pool=list(pool), include_bijection=include_bijection,
                         updated=datetime.datetime.now().isoformat(timespec='seconds'))
    discs = cache['discriminants']
    todo = discriminants_in_range(dmin, dmax)
    done = skipped = failed = 0
    for d in todo:
        key = str(d)
        if (not force) and key in discs and _entry_is_complete(discs[key], include_bijection):
            skipped += 1
            continue
        try:
            discs[key] = compute_disc_entry(d, pool, include_bijection)
            if not discs[key]['success']:
                failed += 1
        except Exception as e:            # never let one d abort the whole run
            discs[key] = {'order': None, 'success': False,
                          'message': f'EXC: {type(e).__name__}: {e}'}
            failed += 1
        done += 1
        if save_every and done % save_every == 0:
            save_cache(cache, path)
            if verbose:
                print(f'  ... {done} computed (at d={d}), {skipped} skipped', flush=True)
    save_cache(cache, path)
    if verbose:
        print(f'done: {done} computed, {skipped} already present, '
              f'{failed} with no rigid l-set; cache holds {len(discs)} discriminants '
              f'-> {path}', flush=True)
    return cache


def _main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        description='Populate / update the rigid l-set + qf-bijection cache.')
    ap.add_argument('--min', type=int, default=-4000,
                    help='most negative discriminant, inclusive (default -4000)')
    ap.add_argument('--max', type=int, default=-3,
                    help='least negative discriminant, inclusive (default -3)')
    ap.add_argument('--path', default=DEFAULT_CACHE, help='cache JSON path')
    ap.add_argument('--no-bijection', action='store_true',
                    help='store only the search data, omit zn_to_qf (smaller file)')
    ap.add_argument('--force', action='store_true',
                    help='recompute entries even if already present')
    ap.add_argument('--save-every', type=int, default=200,
                    help='checkpoint the JSON every N computed entries')
    args = ap.parse_args(argv)
    populate(args.min, args.max, path=args.path,
             include_bijection=not args.no_bijection,
             force=args.force, save_every=args.save_every)


if __name__ == '__main__':
    _main()
