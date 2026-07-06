"""Dry census of SS rigid l-sets for all primes in (1024, 8192].

For each prime, query ss_rigid_lset_info on the disc the assembler uses
(A: -4p, B: -p, C: -4p) and record: case, class-group order, the chosen
ls_rig, descriptor kinds, per-generator (base l, eigenline degree) for every
full-graph generator, the sum descriptor, and the split-l modular-polynomial
coverage for validation.  Dry = no Velu isogenies are computed (~3 min).

    python census.py                 # walk cost (min eigenline degree)
    python census.py --maxcost      # legacy both-direction cost (max degree)

Outputs census_walkcost.json / census_maxcost.json next to this script.
The committed JSONs were produced 2026-07-05; re-run after search changes.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'pycode'))
from ecqf_bij import (ss_rigid_lset_info, _velu_cost, _is_pin, _is_sib,
                      _is_lift, _is_powkey, _desc_base, ssprimes, _modpoly_cache)


def sieve(n):
    isp = bytearray([1]) * (n + 1)
    isp[0] = isp[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if isp[i]:
            isp[i * i::i] = bytearray(len(isp[i * i::i]))
    return [i for i in range(2, n + 1) if isp[i]]


MODPOLY_LS = sorted(set(ssprimes) | {int(k) for k in _modpoly_cache})


def is_split(p, l):
    r = (-p) % l
    return r != 0 and pow(r, (l - 1) // 2, l) == 1


def desc_kind(desc):
    if _is_pin(desc): return 'pin'
    if _is_sib(desc): return 'sib'
    if _is_lift(desc): return 'lift'
    if _is_powkey(desc): return 'powkey'
    if isinstance(desc, tuple): return 'power'
    return 'plain'


def run(walk=True, pmin=1024, pmax=8192):
    def gen_cost(p, desc):
        if _is_sib(desc) or _is_lift(desc):
            return None                     # rational 2-isogeny data, no velu
        l = _desc_base(desc)
        if l == 2:
            return None
        return (l, _velu_cost(0, p, l, walk=walk))

    rows, t0 = [], time.time()
    primes = [p for p in sieve(pmax) if pmin < p]
    for i, p in enumerate(primes):
        d = -p if p % 8 == 7 else -4 * p
        row = {'p': p, 'case': 'A' if p % 4 == 1 else ('B' if p % 8 == 7 else 'C'),
               'd': d}
        info = ss_rigid_lset_info(p, d, walk=walk)
        ls_rig, l_sum = info['ls_rig'], info.get('l_sum')
        gens = [x for x in ls_rig if x != l_sum]
        costs = [c for c in (gen_cost(p, x) for x in gens) if c is not None]
        row.update(order=info['order'], ls_rig=repr(ls_rig),
                   kinds=sorted({desc_kind(x) for x in ls_rig}),
                   success=bool(info.get('success')), gen_costs=costs,
                   max_k=max((k for _, k in costs), default=0),
                   proxy=info['order'] * sum(l * k * k for l, k in costs),
                   ss_ok=not any(desc_kind(x) in ('lift', 'powkey') for x in ls_rig))
        if l_sum is not None:
            row['sum_kind'] = desc_kind(l_sum)
        used = set()
        for x in ls_rig:
            if _is_pin(x):
                if not (_is_sib(x[1]) or _is_lift(x[1])):
                    used.add(_desc_base(x[1]))
                used.update(_desc_base(pk) for pk, c in x[2])
            elif not (_is_sib(x) or _is_lift(x)):
                used.add(_desc_base(x))
        split_ls = [l for l in MODPOLY_LS if l != 2 and is_split(p, l)]
        row['n_split_modpoly'] = len(split_ls)
        row['n_split_unused'] = len([l for l in split_ls if l not in used])
        rows.append(row)
        if (i + 1) % 100 == 0:
            print(f'{i + 1}/{len(primes)}  ({time.time() - t0:.0f}s)', flush=True)
    return rows


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--maxcost', action='store_true',
                    help='legacy both-direction cost model (max degree)')
    args = ap.parse_args()
    walk = not args.maxcost
    rows = run(walk=walk)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'census_walkcost.json' if walk else 'census_maxcost.json')
    json.dump(rows, open(out, 'w'))
    print(f'done: {len(rows)} primes -> {out}')
