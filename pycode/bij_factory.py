"""Batch extension of the ordinary bijection precomputes beyond p < 1024.

Driven by notebooks/bij_factory.ipynb (run it BEFORE phi_factory.ipynb).  The
plan, per the session notes:

1. aps_by_disc(pmin, N): all classes (a, p) with 0 < a < 2 sqrt(p) and p prime
   in (pmin, N], grouped by the exact discriminant d = a^2 - 4p.
2. partition_discs: for each d, look for a rigid spanning set
   (disc_rigid_lset_search, persisted in rigid_lset_cache.json).  Discs split
   into 'spanning' (computable now), 'open' (no rigid l-set from the 15 Atkin
   primes), and 'conductor' (a non-Atkin prime divides cond(d), so the j-side
   vertical scan has no modular polynomial -- a distinct obstruction).
3. compute_and_save: for every (a, p) over a spanning disc, compute the
   j -> form bijection and append it to data/ecqf_ord_pcbij_ext.json, in the
   same format as the shipped ecqf_ord_pcbij_4_1024.json.  Checkpointed per
   prime; re-running resumes.
4. vote_for_new_modpolys: each blocked disc votes for the candidate primes l
   whose (future) modular polynomial would unblock it -- for open discs, the l
   that completes a rigid spanning set; for conductor discs, the unavailable
   conductor primes themselves.  Use the tally to prioritize phi_factory.
5. bootstrap_refresh: the loop tick.  After phi_factory lands new Phi_l's in
   classical_modpolys.json, one call retries the stale failed searches with
   the enlarged pool and mints every newly computable bijection.
"""

import json
import math
import time
from pathlib import Path

from nt import primesBetween, discfac, primefact
from ecqf_bij import disc_rigid_lset_search, ssprimes
from modularpolynomials import modular_prime_pool
from hilbert_crt import modular_set
import rigid_cache

_DATA_DIR = Path(__file__).parent / 'data'
_EXT_BIJ_PATH = _DATA_DIR / 'ecqf_ord_pcbij_ext.json'


#####################################
# Steps 1-2: classes and partition  #
#####################################

def aps_in_range(pmin: int = 1024, pmax: int = 2048) -> list[tuple]:
    """All ordinary classes (a, p): 0 < a < 2 sqrt(p), p prime in (pmin, pmax]."""
    out = []
    for p in primesBetween(pmin + 1, pmax + 1):
        a = 1
        while a * a < 4 * p:
            out.append((a, p))
            a += 1
    return out


def aps_by_disc(pmin: int = 1024, pmax: int = 2048) -> dict:
    """{d: [(a, p), ...]} with d = a^2 - 4p exactly (no closure)."""
    out = {}
    for a, p in aps_in_range(pmin, pmax):
        out.setdefault(a * a - 4 * p, []).append((a, p))
    return out


def _conductor_blockers(d: int) -> list[int]:
    """Primes dividing cond(d) with no modular polynomial available: the
    vertical j-side scan needs one (Atkin format or a classical Phi_l) for
    each conductor prime.

    Exception (the trivial-volcano shortcut in get_ancestor_data_ord): when
    the conductor is EXACTLY one such prime l and the fundamental part has
    class number 1 (d = d0 * l^2, d0 in heeg_js), the l-volcano has a single
    identifiable crater vertex and the vertical scan needs no Phi_l at all."""
    from modularpolynomials import heeg_js
    d0, c = discfac(d)
    mset = modular_set()
    bad = [l for l in primefact(c) if l not in mset]
    if len(bad) == 1 and c == bad[0] and d0 in heeg_js:
        return []
    return bad


def partition_discs(ds: list[int], save_every: int = 200, verbose: bool = True) -> dict:
    """Split ds into {'spanning': [...], 'open': [...], 'conductor': [...]}.

    Runs (and persists) the rigid l-set search for every disc not already in
    rigid_lset_cache.json, so this is the slow, once-per-range step; re-running
    only searches what is missing."""
    cache = rigid_cache.load_cache()
    discs = cache['discriminants']
    out = {'spanning': [], 'open': [], 'conductor': []}
    new = 0
    t0 = time.time()
    for n, d in enumerate(sorted(ds, reverse=True)):
        blockers = _conductor_blockers(d)
        if blockers:
            out['conductor'].append(d)
            continue
        key = str(d)
        if key not in discs:
            try:
                discs[key] = rigid_cache.compute_disc_entry(d)
            except Exception as e:
                discs[key] = {'order': None, 'success': False,
                              'message': f'EXC: {type(e).__name__}: {e}'}
            new += 1
            if new % save_every == 0:
                rigid_cache.save_cache(cache)
                if verbose:
                    print(f'  ... {new} searches done ({n+1}/{len(ds)} discs, '
                          f'{time.time()-t0:.0f}s)', flush=True)
        out['spanning' if discs[key].get('success') else 'open'].append(d)
    if new:
        rigid_cache.save_cache(cache)
    if verbose:
        print(f'partition: {len(out["spanning"])} spanning, {len(out["open"])} open, '
              f'{len(out["conductor"])} conductor-blocked '
              f'({new} new searches, {time.time()-t0:.0f}s)', flush=True)
    return out


def _search_stamp() -> list:
    """Fingerprint of the current search capability: the modular prime pool
    plus feature tags ('sib23' = free sibling generators for conductor-2/3
    kernels of order 2/3; 'lift23' = lifted x_q generators for kernel-1
    descents).  Failed cache entries stamped with the current fingerprint are
    not retried by refresh_blocked_discs."""
    return modular_prime_pool() + ['sib23', 'lift23']


def refresh_blocked_discs(dmin: int = None, dmax: int = -3, save_every: int = 25,
                          verbose: bool = True) -> dict:
    """Retry the cached FAILED rigid searches with the current prime pool.

    partition_discs never re-searches a disc already in rigid_lset_cache.json,
    so discs recorded as open under an older (smaller) pool stay open even
    after new Phi_l's land in classical_modpolys.json.  This recomputes every
    failed cached entry in [dmin, dmax] (default: all of them) with the
    enlarged pool and updates the cache in place.  Conductor-blocked discs are
    not in the cache at all -- they unblock automatically once
    _conductor_blockers stops flagging them.

    Each retried failure is stamped with the search fingerprint it was tried
    under (the prime pool plus search-feature tags), so re-running skips
    entries the current search has already been tried on -- a re-run after
    nothing new has landed is free.  Returns
    {'fixed': [...], 'still_open': [...]}."""
    pool = _search_stamp()
    cache = rigid_cache.load_cache()
    discs = cache['discriminants']
    ds = [int(k) for k, v in discs.items()
          if not v.get('success') and v.get('pool') != pool]
    ds = [d for d in ds if d <= dmax and (dmin is None or d >= dmin)
          and not _conductor_blockers(d)]
    fixed, still = [], []
    t0 = time.time()
    for n, d in enumerate(sorted(ds, reverse=True)):
        try:
            entry = rigid_cache.compute_disc_entry(d)
        except Exception as e:
            entry = {'order': None, 'success': False,
                     'message': f'EXC: {type(e).__name__}: {e}'}
        if not entry.get('success'):
            entry['pool'] = pool
        discs[str(d)] = entry
        (fixed if entry.get('success') else still).append(d)
        if (n + 1) % save_every == 0:
            rigid_cache.save_cache(cache)
            if verbose:
                print(f'  ... {n+1}/{len(ds)} retried, {len(fixed)} unblocked '
                      f'({time.time()-t0:.0f}s)', flush=True)
    rigid_cache.save_cache(cache)
    if verbose:
        print(f'refresh: {len(fixed)} of {len(ds)} unblocked with pool '
              f'{pool} ({time.time()-t0:.0f}s)', flush=True)
    return {'fixed': fixed, 'still_open': still}


def bootstrap_refresh(pmin: int = 1024, pmax: int = 2048, verbose: bool = True) -> dict:
    """One tick of the bootstrap loop: absorb newly minted Phi_l's.

    Run this any time new modular polynomials have landed in
    classical_modpolys.json (phi_factory registers them there itself when a
    target certifies -- there is nothing to register by hand).  It

    1. retries every cached failed rigid search whose last attempt used a
       smaller pool (refresh_blocked_discs), and
    2. re-partitions the range and mints every newly computable bijection
       into the extension store (compute_and_save; this also runs fresh
       searches for discs that just stopped being conductor-blocked).

    Idempotent, and close to free when nothing new has landed.  Returns the
    compute_and_save stats plus 'refresh_fixed' / 'refresh_still_open'."""
    if verbose:
        print(f'bootstrap_refresh with pool {modular_prime_pool()}', flush=True)
    ref = refresh_blocked_discs(verbose=verbose)
    stats = compute_and_save(pmin, pmax, verbose=verbose)
    stats['refresh_fixed'] = ref['fixed']
    stats['refresh_still_open'] = ref['still_open']
    return stats


#########################################
# Step 3: compute and save bijections   #
#########################################

def load_ext_bijections(path=_EXT_BIJ_PATH) -> dict:
    """{(a, p): {j: (A, B, C)}} from the extension store (same format as the
    shipped ecqf_ord_pcbij_4_1024.json)."""
    if not Path(path).exists():
        return {}
    with open(path) as f:
        raw = json.load(f)
    def tup(s):
        return tuple(int(x) for x in s[1:-1].split(','))
    return {tup(k): {int(j): tuple(qf) for j, qf in v.items()} for k, v in raw.items()}


def _save_ext_bijections(store: dict, path=_EXT_BIJ_PATH):
    tmp = str(path) + '.tmp'
    with open(tmp, 'w') as f:
        json.dump({str(ap): {str(j): list(qf) for j, qf in bij.items()}
                   for ap, bij in store.items()}, f)
    Path(tmp).replace(path)


def compute_and_save(pmin: int = 1024, pmax: int = 2048, path=_EXT_BIJ_PATH,
                     verbose: bool = True) -> dict:
    """Compute the bijection for every class over a spanning disc in the range
    and append it to the extension store.  Checkpoints after every prime;
    re-running skips what is already stored.  Returns run statistics."""
    from rigid_cache import ecqf_ord_bij_cached
    by_d = aps_by_disc(pmin, pmax)
    parts = partition_discs(list(by_d), verbose=verbose)
    todo_ds = set(parts['spanning'])
    store = load_ext_bijections(path)
    raw = ({str(ap): {str(j): list(qf) for j, qf in bij.items()}
            for ap, bij in store.items()} if store else {})
    by_p = {}
    for d in todo_ds:
        for a, p in by_d[d]:
            by_p.setdefault(p, []).append((a, p))
    stats = {'computed': 0, 'skipped_existing': 0, 'failed': 0,
             'classes_blocked': sum(len(by_d[d]) for d in parts['open'] + parts['conductor'])}
    t0 = time.time()
    for n, p in enumerate(sorted(by_p)):
        t1 = time.time()
        fresh = 0
        for ap in sorted(by_p[p]):
            if str(ap) in raw:
                stats['skipped_existing'] += 1
                continue
            try:
                bij = ecqf_ord_bij_cached(ap)
            except (ValueError, KeyError) as e:
                stats['failed'] += 1
                if verbose:
                    print(f'  FAILED {ap}: {e}', flush=True)
                continue
            raw[str(ap)] = {str(j): list(qf) for j, qf in bij.items()}
            fresh += 1
        if fresh:
            tmp = str(path) + '.tmp'
            with open(tmp, 'w') as f:
                json.dump(raw, f)
            Path(tmp).replace(path)
            stats['computed'] += fresh
        if verbose:
            done = stats['computed'] + stats['skipped_existing']
            total = sum(len(v) for v in by_p.values())
            print(f'p={p} ({n+1}/{len(by_p)}): {fresh} classes in {time.time()-t1:.1f}s '
                  f'-- {done}/{total} done, {(time.time()-t0)/60:.1f} min elapsed', flush=True)
    stats['elapsed_s'] = time.time() - t0
    return stats


#########################################
# The Hilbert factory                   #
#########################################

def hilbert_factory(save_every: int = 100, verbose: bool = True) -> dict:
    """Compute every Hilbert polynomial the current data certifies.

    Harvests {d: {p: roots of H_d mod p}} from ALL available bijections -- the
    shipped p < 1024 data plus everything in the bij_factory extension store --
    merged with the hilbert_roots_ext cache, then runs the certified-CRT
    reconstruction for every disc whose primes clear its coefficient bound and
    is not already in the library.  Results merge into data/hilbpolys_crt.json
    (checkpointed every save_every discs; re-running skips what is done)."""
    import math
    import hilbert_crt as hc
    from ecqf_tools import ecqf_ord_1K_pc
    t0 = time.time()
    source = {ap: bij for ap, bij in ecqf_ord_1K_pc.items() if ap[0] > 0}
    source.update(load_ext_bijections())
    data = hc.cm_js_by_disc(source)
    lib = hc.hilbert_library()
    if verbose:
        print(f'harvested {len(data)} discs from {len(source)} classes '
              f'({time.time()-t0:.0f}s); library already holds {len(lib)}', flush=True)
    todo = []
    for n, d in enumerate(sorted(data, reverse=True)):
        if d in lib:
            continue
        if sum(math.log2(p) for p in data[d]) > hc.hilbert_bound_log2(d) + 1:
            todo.append(d)
        if verbose and (n + 1) % 4000 == 0:
            print(f'  ... certification scan {n+1}/{len(data)} ({time.time()-t0:.0f}s)',
                  flush=True)
    if verbose:
        print(f'{len(todo)} new discs certifiable; computing '
              f'(checkpoint every {save_every})', flush=True)
    crt_path = hc._DATA_DIR / 'hilbpolys_crt.json'
    out = {}
    if crt_path.exists():
        with open(crt_path) as f:
            out = {int(k): v for k, v in json.load(f).items()}
    stats = {'computed': 0, 'failed': 0}
    for n, d in enumerate(todo):
        try:
            rec = hc.hilbert_via_crt(d, data)
            assert rec['certified'], 'bound not met after selection'
            out[d] = rec['coefs']
            stats['computed'] += 1
        except (ValueError, AssertionError) as e:
            stats['failed'] += 1
            if verbose:
                print(f'  FAILED d={d}: {e}', flush=True)
        if (n + 1) % save_every == 0 or n + 1 == len(todo):
            hc.save_hilbert_library(out)
            if verbose:
                print(f'  {n+1}/{len(todo)} done, library at {len(out)} discs '
                      f'({(time.time()-t0)/60:.1f} min)', flush=True)
    stats['library_size'] = len(out)
    stats['elapsed_s'] = time.time() - t0
    return stats


#########################################
# Step 4: the vote                      #
#########################################

def vote_for_new_modpolys(parts: dict, candidates: tuple = (61, 67, 73, 79, 83, 89, 97, 101),
                          verbose: bool = True) -> dict:
    """Which new modular polynomials Phi_l unblock the most blocked discs?

    Open discs vote for each candidate l that, added to the CURRENT pool (Atkin
    + cached classical Phi_l), yields a rigid spanning set.  Conductor-blocked
    discs vote for their unavailable conductor primes (a Phi_l also fixes the
    vertical scan, via the X1 dispatch).  Returns
    {'votes': {l: [discs]}, 'still_stuck': [discs no single candidate fixes]}."""
    votes = {l: [] for l in candidates}
    stuck = []
    t0 = time.time()
    pool = modular_prime_pool()
    for n, d in enumerate(parts['open']):
        helped = False
        for l in candidates:
            res = disc_rigid_lset_search(d, pool + [l])
            if res['success']:
                votes[l].append(d)
                helped = True
        if not helped:
            stuck.append(d)
        if verbose and (n + 1) % 25 == 0:
            print(f'  ... {n+1}/{len(parts["open"])} open discs polled '
                  f'({time.time()-t0:.0f}s)', flush=True)
    for d in parts['conductor']:
        for l in _conductor_blockers(d):
            if l in votes:
                votes[l].append(d)
    if verbose:
        for l in sorted(votes, key=lambda l: -len(votes[l])):
            if votes[l]:
                print(f'Phi_{l}: unblocks {len(votes[l])} discs')
        print(f'{len(stuck)} open discs need more than one new prime')
    return {'votes': votes, 'still_stuck': stuck}
