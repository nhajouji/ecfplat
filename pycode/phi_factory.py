"""Batch production of classical modular polynomials Phi_p -- the bootstrap driver.

Designed for a long unattended run on a separate machine (see
notebooks/phi_factory.ipynb): for each prime ell above the precomputed range it
extends the rigid-l-set cache, computes the j <-> form bijections for every
class at ell (skipping the rare rigid gaps), extracts p-isogenous pairs for
every unfinished target p, solves for Phi_p mod ell, and folds the solution
into a running CRT.  A target finishes when its modulus provably exceeds the
Broker-Sutherland height bound; the finished matrix must pass the Kronecker
congruence and is then registered into classical_modpolys.json.

Everything checkpoints: the rigid cache persists per disc, and the CRT state
(data/phi_factory_state.json) persists after every prime, so the run can be
interrupted and resumed with the same 'run all'.

Why this needs large ell: a full solve mod ell takes ~(p+1)(p+2)/2 independent
equations and a prime supplies ~0.55*ell pairs, so p = 37..67 only start
collecting primes at ell ~ 2500..4500 -- beyond the shipped ell < 1024 data.
No special-value relations are required (they lower the usable-ell threshold
when present in the Hilbert library, and are picked up automatically).
"""

import json
import math
import time
from pathlib import Path

from nt import primeQ, primesBetween, crt_pair
from qfs import qf_isogs_hor
from ecqf_tools import ecqf_ord_1K_pc
from modularpolynomials import classical_modpoly, register_modpoly, _modpoly_cache
from modpoly_crt import (phi_monomials, phi_diagonal, special_value_rows,
                         solve_phi_from_pairs, isogenous_pairs_mod_ell,
                         modpoly_height_bound_log2, kronecker_check)
from hilbert_crt import hilbert_library
import rigid_cache

_DATA_DIR = Path(__file__).parent / 'data'
_STATE_PATH = _DATA_DIR / 'phi_factory_state.json'
_SHIPPED_ELLS = frozenset(l0 for a, l0 in ecqf_ord_1K_pc if a > 0)

WANTED = (37, 43, 53, 61, 67)                 # the non-Atkin primes below 72


def factory_targets(lmax_p: int = 71) -> list[int]:
    """Primes p <= lmax_p with no classical modular polynomial in the cache."""
    return [p for p in primesBetween(1, lmax_p + 1) if p not in _modpoly_cache]


##############################
# Checkpointed CRT state     #
##############################

def load_state(path=_STATE_PATH) -> dict:
    if Path(path).exists():
        with open(path) as f:
            raw = json.load(f)
        return {int(p): st for p, st in raw.items()}
    return {}


def save_state(state: dict, path=_STATE_PATH):
    tmp = str(path) + '.tmp'
    with open(tmp, 'w') as f:
        json.dump({str(p): st for p, st in state.items()}, f)
    Path(tmp).replace(path)


def _fold_solution(st: dict, sol: dict, ell: int, mons: list):
    """CRT-accumulate one mod-ell solution into the target's running state."""
    vec = [sol[m] for m in mons]
    if st['modulus'] == 1:
        st['crt'] = [v if 2 * v <= ell else v - ell for v in vec]
        st['modulus'] = ell
    else:
        M = st['modulus']
        st['crt'] = [crt_pair((c, M), (v, ell))[0] for c, v in zip(st['crt'], vec)]
        st['modulus'] = M * ell
    st['ells'].append(ell)


##############################
# Pairs at a fresh prime     #
##############################

_POP_FLOOR = [0]                               # most negative disc the rigid cache covers

def _ensure_rigid(ell: int, headroom: int = 8192):
    """Extend the rigid-l-set cache to cover every class disc at ell, in chunks:
    populating past the immediate need (headroom) keeps the big cache file from
    being rewritten on every prime."""
    need = -4 * ell
    if _POP_FLOOR[0] <= need:
        return
    target = need - headroom
    print(f'  [rigid cache: extending to {target} ...]', flush=True)
    rigid_cache.populate(target, save_every=500, verbose=True)
    _POP_FLOOR[0] = target


_EXT_BY_P = None                               # bij_factory store, grouped by prime

def _ext_bijections_at(ell: int):
    global _EXT_BY_P
    if _EXT_BY_P is None:
        from bij_factory import load_ext_bijections
        _EXT_BY_P = {}
        for (a, p), bij in load_ext_bijections().items():
            _EXT_BY_P.setdefault(p, []).append(bij)
    return _EXT_BY_P.get(ell)


def class_bijections_at_ell(ell: int, verbose: bool = False,
                            headroom: int = 8192) -> list[dict]:
    """The j -> form bijections of every ordinary class at ell, skipping classes
    whose disc has no rigid l-set (the gap discs).  Uses the bij_factory
    extension store when it covers ell (run notebooks/bij_factory.ipynb first
    and this costs nothing); otherwise computes on the fly, extending the rigid
    cache on disk as a side effect."""
    stored = _ext_bijections_at(ell)
    if stored:
        return stored
    _ensure_rigid(ell, headroom)
    from rigid_cache import ecqf_ord_bij_cached
    out = []
    a = 1
    while a * a < 4 * ell:
        if (a * a - 4 * ell) % 4 in (0, 1):
            try:
                out.append(ecqf_ord_bij_cached((a, ell)))
            except (ValueError, KeyError) as e:
                if verbose:
                    print(f'    skip class ({a}, {ell}): {e}')
        a += 1
    return out


def pairs_from_bijections(bijections: list[dict], p: int) -> list[tuple]:
    """p-isogenous (j1, j2) pairs over F_ell from per-class bijections."""
    pairs = set()
    for bij in bijections:
        qf_to_j = {qf: j for j, qf in bij.items()}
        for j1, qf in bij.items():
            for qf2 in qf_isogs_hor(qf, p):
                j2 = qf_to_j.get(tuple(qf2))
                if j2 is not None:
                    pairs.add((min(j1, j2), max(j1, j2)))
    return sorted(pairs)


##############################
# The driver                 #
##############################

def _finalize(p: int, st: dict, mons: list, save: bool = True) -> bool:
    """Build the matrix from the CRT state, Kronecker-check, register + persist."""
    M = [[0] * (p + 2) for _ in range(p + 2)]
    M[p + 1][0] = M[0][p + 1] = 1
    for (i, j), c in zip(mons, st['crt']):
        M[i][j] = M[j][i] = c
    if not kronecker_check(p, M):
        return False
    register_modpoly(p, M, save=save)
    return True


def run_factory(targets: list[int] = None, lmin: int = 1024, lmax: int = 16384,
                pair_cap: float = 1.5, save: bool = True, use_1k_data: bool = True,
                headroom: int = 8192):
    """Compute Phi_p for every target, checkpointing after each prime.

    Prints one status line per prime with per-target progress (bits collected vs
    the certification bound).  Re-running resumes from the checkpoint; finished
    targets are skipped.  pair_cap limits the pair equations per system to
    pair_cap * #unknowns (solver speed)."""
    if targets is None:
        targets = factory_targets()
    state = load_state()
    hilb = hilbert_library()
    work = {}
    for p in targets:
        if p in _modpoly_cache:
            continue
        mons = phi_monomials(p)
        st = state.setdefault(p, {'crt': [], 'modulus': 1, 'ells': [], 'done': False})
        if st['done']:
            continue
        extra = special_value_rows(p, hilb)
        work[p] = {'mons': mons, 'diag': phi_diagonal(p, hilb).int_coefs(),
                   'extra': extra if extra['rows'] else None,
                   'bound': modpoly_height_bound_log2(p),
                   'min_rows': len(mons) - (2 * p + 1) - len(extra['rows']),
                   'st': st}
        print(f'target p={p}: {len(mons)} unknowns, bound {work[p]["bound"]:.0f} bits, '
              f'{math.log2(st["modulus"]):.0f} bits done, '
              f'{len(extra["used"])} special-value families', flush=True)
    if not work:
        print('nothing to do: all targets already in classical_modpolys.json')
        return

    ells = [l for l in primesBetween(lmin, lmax + 1)]
    if use_1k_data:                            # the shipped ell < 1024 bijections are free
        ells = sorted({l0 for a, l0 in ecqf_ord_1K_pc if a > 0} | set(ells))
    for ell in ells:
        todo = {p: w for p, w in work.items() if not w['st']['done']
                and ell not in w['st']['ells'] and ell != p}
        if not todo:
            continue                           # nothing to do at THIS prime only
        fresh = ell not in _SHIPPED_ELLS
        t0 = time.time()
        bijections = class_bijections_at_ell(ell, headroom=headroom) if fresh else None
        line = []
        for p, w in sorted(todo.items()):
            st = w['st']
            pairs = (pairs_from_bijections(bijections, p) if fresh
                     else isogenous_pairs_mod_ell(p, ell))
            if len(pairs) < w['min_rows']:     # cannot reach full rank: skip cheaply
                continue
            cap = int(pair_cap * len(w['mons']))
            sol = solve_phi_from_pairs(p, ell, w['diag'], pairs[:cap], w['extra'])
            if sol is None:
                continue
            _fold_solution(st, sol, ell, w['mons'])
            bits = math.log2(st['modulus'])
            if bits > w['bound'] + 1:
                st['done'] = _finalize(p, st, w['mons'], save=save)
                line.append(f'p={p} {"FINISHED" if st["done"] else "KRONECKER FAIL"}')
            else:
                line.append(f'p={p} {bits:.0f}/{w["bound"]:.0f}')
        if save:
            save_state(state)
        if line:
            print(f'ell={ell} ({time.time()-t0:.0f}s): ' + '  '.join(line), flush=True)
        if all(w['st']['done'] for w in work.values()):
            break

    print('\nsummary:', flush=True)
    for p, w in sorted(work.items()):
        st = w['st']
        if st['done']:
            msg = 'DONE, saved to classical_modpolys.json'
        else:
            msg = (f'in progress -- {math.log2(st["modulus"]):.0f} / {w["bound"]:.0f} bits '
                   f'from {len(st["ells"])} primes (rerun with larger lmax to continue)')
        print(f'  p={p}: {msg}', flush=True)
