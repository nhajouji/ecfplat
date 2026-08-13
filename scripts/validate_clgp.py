"""Validation sweep for pycode/clgp.py (class groups via Gauss composition).

Synchronous and small by design (~1 min).  Six batteries:
  A. group axioms + |Cl(d)| == clgr_size_gen, all discriminants -4 >= d >= -800
  B. composition vs isogeny cycles: <x_l> == qf_isog_cycle pointwise, d >= -400
  C. 2-torsion count vs genus theory (clgr2_size), fundamental d
  D. basis / invariant factors (+ known spot checks)
  E. prime buckets really represent their primes (brute-force evaluation)
  F. rigid_lset contract + mutual check against disc_rigid_lset_search

Run:  python scripts/validate_clgp.py
"""

import random
import sys
import time
from math import lcm
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'pycode'))

from nt import discfac, primesBetween
from identities import clgr_size_gen, clgr2_size
from alg_classes import Subgroup, element_order
from qfs import (class_group_id, class_group_inv, get_qfs_strict,
                 qf_mod_gamma, qf_isog_cycle, qf_evs_inrange)
from clgp import (class_group, qf_compose, class_from_prime, primes_by_class,
                  rigid_lset)

D_MIN_AXIOMS = -800
D_MIN_CYCLES = -400
D_MIN_PRIMES = -200

FAILURES = []

def check(cond, label):
    if not cond:
        FAILURES.append(label)
        print(f'  FAIL: {label}')

def discs(dmin):
    return [d for d in range(-3, dmin - 1, -1) if d % 4 in (0, 1)]


# --- A: axioms and order ---
t0 = time.time()
rng = random.Random(0)
n_comp = 0
for d in discs(D_MIN_AXIOMS):
    G = class_group(d)
    elems = get_qfs_strict(d)
    check(len(elems) == clgr_size_gen(d) == G.order, f'A |Cl({d})|')
    check(G.zero_element in elems, f'A identity in elements d={d}')
    for x in elems:
        check(G.add_elements(x, class_group_inv(x)) == G.zero_element,
              f'A inverse d={d} x={x}')
        check(G.add_elements(x, G.zero_element) == x, f'A identity d={d} x={x}')
    eset = set(elems)
    for x in elems:
        for y in elems:
            s = G.add_elements(x, y)
            n_comp += 1
            check(s in eset, f'A closure d={d} {x}+{y}')
            check(s == G.add_elements(y, x), f'A commutativity d={d} {x},{y}')
    for _ in range(min(30, len(elems) ** 3)):
        x, y, z = (rng.choice(elems) for _ in range(3))
        check(G.add_elements(G.add_elements(x, y), z)
              == G.add_elements(x, G.add_elements(y, z)),
              f'A associativity d={d} {x},{y},{z}')
print(f'A: axioms on {len(discs(D_MIN_AXIOMS))} discs, {n_comp} compositions '
      f'({time.time()-t0:.1f}s)')

# --- B: composition vs isogeny cycles ---
t0 = time.time()
n_checked = 0
for d in discs(D_MIN_CYCLES):
    G = class_group(d)
    qf0 = qf_mod_gamma(class_group_id(d))
    for l in primesBetween(2, 40):
        x = class_from_prime(d, l)
        if x is None or x == qf0:
            continue
        cyc = qf_isog_cycle(qf0, l)
        o = element_order(G, x, G.order)
        check(len(set(cyc)) == o, f'B cycle length d={d} l={l}')
        check(set(cyc) == Subgroup(G, [x]).elements, f'B cycle set d={d} l={l}')
        # pointwise: after orienting, cyc[k] == k*(+-x)
        xo = x if cyc[1 % len(cyc)] == x else class_group_inv(x)
        walk, ok = qf0, True
        for k in range(len(set(cyc))):
            if cyc[k] != walk:
                ok = False
                break
            walk = G.add_elements(walk, xo)
        check(ok, f'B pointwise walk d={d} l={l}')
        n_checked += 1
print(f'B: {n_checked} (d, l) cycles matched against composition '
      f'({time.time()-t0:.1f}s)')

# --- C: 2-torsion vs genus theory ---
n_checked = 0
for d in discs(D_MIN_AXIOMS):
    if abs(d) < 13 or discfac(d)[1] != 1:
        continue
    G = class_group(d)
    t2 = sum(1 for x in get_qfs_strict(d)
             if G.add_elements(x, x) == G.zero_element)
    check(t2 == clgr2_size(d), f'C |Cl[2]| d={d}: {t2} vs {clgr2_size(d)}')
    n_checked += 1
print(f'C: 2-torsion count vs genus theory on {n_checked} fundamental discs')

# --- D: basis and invariant factors ---
t0 = time.time()
for d in discs(D_MIN_AXIOMS):
    G = class_group(d)
    basis = G.basis()
    prod = 1
    for _, o in basis:
        prod *= o
    check(prod == G.order, f'D basis order product d={d}')
    check(len(Subgroup(G, [x.vec for x, _ in basis])) == G.order,
          f'D basis spans d={d}')
    inv = G.structure()
    check(all(inv[i + 1] % inv[i] == 0 for i in range(len(inv) - 1)),
          f'D divisibility chain d={d}')
    expo = lcm(*(element_order(G, x, G.order) for x in get_qfs_strict(d)))
    check((inv[-1] if inv else 1) == expo, f'D exponent d={d}')
for d, want in [(-39, (4,)), (-47, (5,)), (-56, (4,)), (-84, (2, 2)),
                (-95, (8,)), (-480, (2, 2, 2))]:
    got = class_group(d).structure()
    check(got == want, f'D known structure d={d}: got {got} want {want}')
print(f'D: bases + invariant factors on {len(discs(D_MIN_AXIOMS))} discs '
      f'({time.time()-t0:.1f}s)')

# --- E: prime buckets brute-checked by evaluation ---
n_checked = 0
for d in discs(D_MIN_PRIMES):
    buckets = primes_by_class(d, 200)
    for qf, ls in buckets.items():
        if not ls:
            continue
        l = ls[0]
        check(l in qf_evs_inrange(qf, 20), f'E {qf} represents {l} (d={d})')
        n_checked += 1
print(f'E: {n_checked} smallest-prime buckets brute-checked by evaluation')

# --- F: rigid_lset contract + mutual check vs disc_rigid_lset_search ---
t0 = time.time()
n_ok = 0
for d in discs(D_MIN_AXIOMS):
    G = class_group(d)
    out = rigid_lset(G)
    check(out is not None, f'F rigid_lset succeeds d={d}')
    if out is None:
        continue
    basis = out['basis']
    prod = 1
    for _, _, o in basis:
        prod *= o
    check(prod == G.order, f'F basis product d={d}')
    if basis:
        check(len(Subgroup(G, [qf for _, qf, _ in basis])) == G.order,
              f'F basis spans d={d}')
    for l, qf, o in basis:
        check(class_from_prime(d, l) in (qf, class_group_inv(qf)),
              f'F prime {l} represents its generator d={d}')
        check(element_order(G, qf, G.order) == o, f'F order of {l}-gen d={d}')
    if out['pin'] is not None:
        xpin = class_from_prime(d, out['pin'])
        sums = {G.zero_element}
        for _, qf, _o in (b for b in basis if b[2] > 2):
            sums = ({G.add_elements(s, qf) for s in sums}
                    | {G.add_elements(s, class_group_inv(qf)) for s in sums})
        check(xpin in sums or class_group_inv(xpin) in sums,
              f'F pin is a signed sum d={d}')
    n_ok += 1
t_new = time.time() - t0
print(f'F: rigid_lset contract verified on {n_ok} discs ({t_new:.1f}s)')

# mutual check + timing against the production search on a sample
from ecqf_bij import disc_rigid_lset_search
sample = [-231, -255, -420, -480, -560, -644, -3299, -6272]
t_old = t_grp = 0.0
for d in sample:
    t0 = time.time()
    ref = disc_rigid_lset_search(d)
    t_old += time.time() - t0
    t0 = time.time()
    G = class_group(d)
    mine = rigid_lset(G)
    t_grp += time.time() - t0
    if not ref['success'] or any(not isinstance(x, int) for x in ref['ls_rig']):
        print(f'  note: d={d} production search used descriptors/failed; skipped')
        continue
    # validate the PRODUCTION result with the group law
    for l, o in zip(ref['ls'], ref['ns']):
        x = class_from_prime(d, l)
        check(x is not None and element_order(G, x, G.order) == o,
              f'F ref order of x_{l} d={d}')
    span = Subgroup(G, [class_from_prime(d, l) for l in ref['ls']])
    check(len(span) == G.order, f'F ref basis spans d={d}')
    if ref['l_sum'] is not None:
        xpin = class_from_prime(d, ref['l_sum'])
        sums = {G.zero_element}
        for l in ref['ls_basis']:
            qf = class_from_prime(d, l)
            sums = ({G.add_elements(s, qf) for s in sums}
                    | {G.add_elements(s, class_group_inv(qf)) for s in sums})
        check(xpin in sums or class_group_inv(xpin) in sums,
              f'F ref pin is a signed sum d={d}')
    check(mine is not None, f'F group-law rigid_lset succeeds d={d}')
print(f'F: mutual check on {len(sample)} discs -- production search '
      f'{t_old:.1f}s vs group-law {t_grp:.2f}s')

print()
if FAILURES:
    print(f'{len(FAILURES)} FAILURES')
    sys.exit(1)
print('ALL CHECKS PASSED')
