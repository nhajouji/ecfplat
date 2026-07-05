from nt import primesBetween,discfac,quad_rec
from qfs import *
from ecfp import trfr_to_js, js_to_fg
from graph_tools import compute_bijection_zn, cycle_from_neighbor_data
from modularpolynomials import (eval_atk, small_bij_check, atkin_polys_dict,
                                modpoly_roots_among, modular_prime_pool,
                                _modpoly_cache, heeg_js)
from misctools import *

import itertools
import json
from functools import lru_cache
from math import gcd
from pathlib import Path

ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]

_DATA_DIR = Path(__file__).parent / 'data'

###################################
# Loading precomputed l-set data  #
###################################
# Precomputed disc_rigid_lset_search(d) output keyed by discriminant d, so the
# (moderately expensive) class-group l-set search need not be rerun for cached d.
# Mirrors the precomputed-bijection loading in ecqf_tools.py.  Read it through
# disc_ldata(d) (defined below) for an automatic fall-back to a live search.

def _ldata_to_tuples(entry):
    out = dict(entry)
    for key in ('ls_rig', 'ls_basis', 'ls_2tors', 'ls', 'ns', 'best'):
        v = out.get(key)
        if isinstance(v, list):                     # power descriptors [l,k] -> (l,k)
            out[key] = tuple(tuple(x) if isinstance(x, list) else x for x in v)
    if isinstance(out.get('l_sum'), list):
        out['l_sum'] = tuple(out['l_sum'])
    return out

try:
    with open(_DATA_DIR / 'qf_ldata.json') as f:
        _qf_ldata_loaded = json.load(f)
    qf_ldata = {int(d): _ldata_to_tuples(_qf_ldata_loaded[d]) for d in _qf_ldata_loaded}
except FileNotFoundError:
    qf_ldata = {}     # not generated yet -> disc_ldata falls back to a live search
                      # (regenerate with: python ldata_cache.py --min -4096)


def qf_reps_pm(d:int):
    return [qf for qf in get_qfs_strict(d) if qf[1]>=0]

def gen_qfs_d_ls(d,lset):
    qf0 = class_group_id(d)
    qfs = [qf0]
    for l in lset:
        qfs_new = []
        for qf in qfs:
            qfs_new+=qf_isog_cycle(qf,l)
        qfs = list(set(qfs_new))
    return qfs


                ####################################
                # Rigid l-set search (from scratch) #
                ####################################
#
# Goal: given the class group A of discriminant d and a pool of candidate
# primes `ls` (each prime l represents an element x_l in A via the horizontal
# isogeny action), find a subset whose elements form a *rigid spanning set*,
# packaged so it can be handed straight to `ecqf_full_bijection_ord` /
# `compute_bijection_zn`.
#
# What the downstream cube construction in graph_tools needs from the tuple it
# is given:
#   * The elements must be INDEPENDENT generators, i.e. a direct-product basis
#     A = prod <x_l>.  For a generating set this is equivalent to
#         prod_l ord(x_l) == |A|
#     (a generating set always over-counts; equality <=> direct product).
#     This is strictly stronger than "irredundant" but is what the labelling
#     algorithm requires, and it implies irredundancy.
#   * Generators of order 2 are handled separately (extend_bijection_zn2) and
#     just need to be independent 2-torsion; no orientation data is needed.
#   * If there are m >= 2 generators of order > 2, the cube needs one extra
#     "pinning" element to fix the relative orientation of the m cycles.
#     compute_bijection_zn uses the *last* long prime for this, and the tree
#     search only closes up if that prime's element is a signed sum
#         +- x_1 +- ... +- x_m
#     of the other long generators -- i.e. exactly the "sum of the order > 2
#     generators" x* from the notebook.  With m <= 1 no pinning element is
#     needed.
#
# So a valid rigid tuple is
#     ls_rig = (b_1, ..., b_m, l_sum, c_1, ..., c_k)
# with b_i the independent order > 2 generators, l_sum a prime representing a
# signed sum of them (omitted when m <= 1), and c_j independent order-2
# generators (their position is irrelevant -- l2_split pulls them out by order).


def _l_element_order(d, l, qf0=None):
    """Order of the element x_l = length of its horizontal isogeny cycle."""
    if qf0 is None:
        qf0 = class_group_id(d)
    return len(qf_isog_cycle(qf0, l))


def _rig_candidates(d, ls, qf0):
    """Map each usable prime to ord(x_l), dropping primes with trivial x_l."""
    cand = {}
    for l in ls:
        nontriv = [qf for qf in qf_isogs_hor(qf0, l) if qf != qf0]
        if len(nontriv) == 0:
            continue                       # x_l is the identity -> useless
        o = _l_element_order(d, l, qf0)
        if o > 1:
            cand[l] = o
    return cand


# ----------------------------------------------------------------------------
# Generator descriptors
# ----------------------------------------------------------------------------
# A "direction" is either a bare prime l (a horizontal l-isogeny) or a power
# (l, k), meaning the element k * x_l.  A power's neighbour data is read off the
# l-cycle by stepping k -- no new isogeny data is needed.  Labels stay plain
# ints when k == 1, so prime-only results are byte-identical to before.

def _divisors(n):
    return [k for k in range(1, n + 1) if n % k == 0]

def _desc_base(desc):
    return desc[0] if isinstance(desc, tuple) else desc

# Sibling ("free prime") descriptors ('sib', q), q in {2, 3}.  When q divides
# the conductor, the fibers of the descent Cl(O) -> Cl(O_{c/q}) -- the classes
# sharing a q-parent, "q-siblings" -- are cosets of its kernel.  When that
# kernel has order 2 or 3, the UNORDERED fiber data already matches the kernel
# generator's cycle structure (any arrangement of 2 or 3 elements is a valid
# cycle; with 4 the pairings abcd vs acbd are inequivalent), so the fibers can
# serve as a generator direction in a rigid l-set with NO modular polynomial
# beyond Phi_q itself: on the j side the fibers are read off the q-ancestor
# data (Atkin q = 2, 3), on the qf side off qf_isogs_des.

def _is_sib(desc):
    return isinstance(desc, tuple) and desc[0] == 'sib'

def qf_q_siblings(qf, q):
    """The full fiber (qf included) of classes sharing qf's q-parent."""
    d = qf_disc(qf)
    par = qf_isog_parent(tuple(qf), q)
    return list(dict.fromkeys(x for x in qf_isogs_des(par, q) if qf_disc(x) == d))

def sib_kernel_order(d, q):
    """|ker(Cl(d) -> Cl(d/q^2))| = size of the sibling fibers; 0 when q does
    not divide cond(d)."""
    if discfac(d)[1] % q:
        return 0
    return len(qf_q_siblings(class_group_id(d), q))

def _sib_candidates(d, qf0=None):
    """Usable free-generator candidates [(('sib', q), ('sib', q), order)]."""
    out = []
    for q in (2, 3):
        k = sib_kernel_order(d, q)
        if k in (2, 3):
            out.append((('sib', q), ('sib', q), k))
    return out

# Lifted generators ('lift', q): when the kernel of Cl(d) -> Cl(d/q^2) is
# TRIVIAL (q splits one level up, so every parent class has exactly one child
# and the parent map is a bijection), the crater's x_q-action conjugates down
# to a floor generator: j1 -> j2 iff their q-parents are horizontally
# q-isogenous.  This makes q itself usable as a generator even though q
# divides the conductor -- e.g. d = -4*71 gets its complete bijection from the
# 2-isogeny graph alone.  Equivalent to computing the bijection upstairs and
# extending down the unique descents.

def qf_q_child(qf_par, q, d):
    """The unique disc-d class below qf_par (only valid when the kernel of the
    q-descent is trivial)."""
    kids = list(dict.fromkeys(x for x in qf_isogs_des(qf_par, q) if qf_disc(x) == d))
    assert len(kids) == 1, f'{qf_par}: descent to disc {d} is not kernel-1'
    return kids[0]

def qf_lift_nbrs(qf, q):
    """The +- neighbours of qf under the lifted x_q: the unique children of
    its q-parent's horizontal q-neighbours."""
    d = qf_disc(qf)
    par = qf_isog_parent(tuple(qf), q)
    return list(dict.fromkeys(
        qf_q_child(P, q, d) for P in dict.fromkeys(qf_isogs_hor(par, q))))

def _lift_candidates(d, qf0=None):
    """Usable lifted-generator candidates [(('lift', q), ('lift', q), order)]."""
    if qf0 is None:
        qf0 = class_group_id(d)
    out = []
    for q in (2, 3):
        if sib_kernel_order(d, q) == 1:                 # kernel-1 descent
            n = len(qf_isog_cycle(qf_isog_parent(qf0, q), q))
            if n > 1:
                out.append((('lift', q), ('lift', q), n))
    return out

def _is_lift(desc):
    return isinstance(desc, tuple) and desc[0] == 'lift'

# Relaxed pin descriptors ('pin', pdesc, pairs, taus): the pinning element is
# sum(c_i * b_i) + tau instead of the plain sum of the long generators, with
# full support, 0 < c_i < ord(b_i)/2, and tau a word in the basis order-2
# generators.  pdesc is the candidate representing it; pairs[i] = (key of the
# c_i-step graph of long generator i, c_i); taus = the order-2 descriptors of
# the tau word.  graph_tools.compute_bijection_zn interprets the label; the
# data builders below expand it into the needed neighbour graphs.

def _is_pin(desc):
    return isinstance(desc, tuple) and len(desc) == 4 and desc[0] == 'pin'

def _pow_desc(desc, c):
    """The neighbour-data key for the c-step of a generator direction."""
    if c == 1:
        return desc
    if _is_sib(desc) or _is_lift(desc):
        return ('pow', desc, c)
    if isinstance(desc, tuple):
        return (desc[0], desc[1] * c)
    return (desc, c)

def _is_powkey(desc):
    return isinstance(desc, tuple) and len(desc) == 3 and desc[0] == 'pow'

def _desc_cycle(qf, desc):
    """Cyclic orbit of qf under the descriptor direction (qf side)."""
    if _is_sib(desc):
        return qf_q_siblings(qf, desc[1])
    if _is_lift(desc):
        q = desc[1]
        d = qf_disc(qf)
        return [qf_q_child(P, q, d)
                for P in qf_isog_cycle(qf_isog_parent(tuple(qf), q), q)]
    return qf_isog_cycle_power(qf, desc) if isinstance(desc, tuple) else qf_isog_cycle(qf, desc)

def _desc_order(d, desc, qf0=None):
    """Order of the element a descriptor represents."""
    if qf0 is None:
        qf0 = class_group_id(d)
    if _is_sib(desc):
        return sib_kernel_order(d, desc[1])
    if _is_lift(desc):
        return len(qf_isog_cycle(qf_isog_parent(qf0, desc[1]), desc[1]))
    if isinstance(desc, tuple):
        l, k = desc
        n = len(qf_isog_cycle(qf0, l))
        return n // gcd(n, k)
    return len(qf_isog_cycle(qf0, desc))

def gen_from_descs(d, descs):
    """Subgroup generated by a list of descriptor directions (qf side).

    For a list of bare primes this coincides with gen_qfs_d_ls."""
    qf0 = class_group_id(d)
    qfs = {qf0}
    for desc in descs:
        nxt = set()
        for qf in qfs:
            for q in _desc_cycle(qf, desc):
                nxt.add(q)
        qfs = nxt
    return qfs

def _desc_qf_pm(qf, desc):
    """The de-duplicated +- neighbours of qf under the descriptor (qf side)."""
    if _is_sib(desc):
        return [x for x in qf_q_siblings(qf, desc[1]) if x != qf]
    if _is_lift(desc):
        return qf_lift_nbrs(qf, desc[1])
    if isinstance(desc, tuple):
        l, k = desc
        cyc = qf_isog_cycle(qf, l)
        n = len(cyc)
        return list(dict.fromkeys([cyc[k % n], cyc[(-k) % n]]))
    return list(dict.fromkeys(qf_isogs_hor(qf, desc)))

def _kstep_nbrdata(base_nbr, k):
    """Turn an l-cycle graph {v:[v+-x_l]} into the (l,k) graph {v:[v+-k x_l]}.

    Orientation-free: for each v it returns the unordered pair {v+k x_l, v-k x_l}
    read off v's cycle, so it works identically on the qf and j sides."""
    out = {}
    for v in base_nbr:
        cyc = cycle_from_neighbor_data(v, base_nbr)
        n = len(cyc)
        out[v] = list(dict.fromkeys([cyc[k % n], cyc[(-k) % n]]))
    return out

def _signed_sum_corners(d, long_descs, qf0):
    """All cube far-corners  +- x_1 +- ... +- x_m  for the long generators.

    These are the endpoints the cube tree-search produces; the pinning generator
    must represent one of them."""
    frontier = {qf0}
    for desc in long_descs:
        nxt = set()
        for qf in frontier:
            for qf1 in _desc_qf_pm(qf, desc):
                nxt.add(qf1)
        frontier = nxt
    return frontier

def _iter_bases_desc(d, cands, target):
    """DFS over ALL independent generating sets from descriptor candidates,
    yielded lazily in search order (the first yield reproduces the historical
    single-basis search).

    cands is a list of (desc, base_prime, order) sorted by order descending; at
    most one descriptor per base prime is used (all powers of x_l share <x_l>).
    Yielding every basis lets the rigidity step backtrack: a basis whose sum
    element no candidate represents is not the end -- another basis of the
    same group may admit a pin."""
    def dfs(start, S, used, g):
        if g == target:
            yield list(S)
            return
        for i in range(start, len(cands)):
            desc, base, m = cands[i]
            if base in used:
                continue
            if len(gen_from_descs(d, [c[0] for c in S] + [desc])) == g * m:
                yield from dfs(i + 1, S + [cands[i]], used | {base}, g * m)
    yield from dfs(0, [], set(), 1)

def _best_independent_partial(d, cand_sorted, cand):
    """Largest independent prime-only set -- used to report a near miss."""
    best = []
    best_g = 1
    def dfs(start, S, g):
        nonlocal best, best_g
        if g > best_g:
            best_g, best = g, list(S)
        for i in range(start, len(cand_sorted)):
            l = cand_sorted[i]
            g2 = len(gen_qfs_d_ls(d, S + [l]))
            if g2 == g * cand[l]:
                dfs(i + 1, S + [l], g2)
    dfs(0, [], 1)
    return best

def _prime_candidates(cand):
    """Bare-prime descriptors, highest order first (reproduces the old search)."""
    return [(l, l, cand[l]) for l in sorted(cand, key=lambda l: cand[l], reverse=True)]

def _power_candidates(cand):
    """Bare primes plus their proper-power descriptors (l, k), highest order first."""
    cands = []
    for l in cand:
        n = cand[l]
        cands.append((l, l, n))                               # the prime itself
        for m in sorted((mm for mm in _divisors(n) if 1 < mm < n), reverse=True):
            cands.append(((l, n // m), l, m))                 # element of order m
    cands.sort(key=lambda c: c[2], reverse=True)
    return cands

def _rigid_from_candidates(d, qf0, order, cands):
    """Find an independent basis among cands and assemble the rigid l-set.

    Returns ('ok', longs, twos, l_sum, ls_rig), ('rigidfail', longs, twos), or
    None (no independent basis among these candidates).  longs/twos are lists of
    (desc, base, order) triples.

    Backtracks over bases: when a basis needs a pinning sum element that no
    candidate represents, the DFS continues to the next basis -- different
    bases of the same group have different sum elements, and one of them may
    be represented (e.g. d = -6272: (3, 11) has no pin in the pool, but
    (3, 43) is pinned by 17).  The first rigidfail is kept for reporting."""
    first_fail = None
    for basis in _iter_bases_desc(d, cands, order):
        longs = sorted([c for c in basis if c[2] > 2], key=lambda c: c[2], reverse=True)
        twos = [c for c in basis if c[2] == 2]
        l_sum = None
        if len(longs) >= 2:
            corners = _signed_sum_corners(d, [c[0] for c in longs], qf0)
            used = {c[1] for c in basis}
            for desc, base, m in cands:
                if base in used:
                    continue
                if any(qf in corners for qf in _desc_qf_pm(qf0, desc) if qf != qf0):
                    l_sum = desc
                    break
            if l_sum is None:
                if first_fail is None:
                    first_fail = ('rigidfail', longs, twos)
                continue
        ls_rig = (tuple(c[0] for c in longs)
                  + ((l_sum,) if l_sum is not None else ())
                  + tuple(c[0] for c in twos))
        return ('ok', longs, twos, l_sum, ls_rig)
    return first_fail


def _rigid_relaxed_pin(d, qf0, order, cands):
    """The relaxed rigidity pass: the pinning element may be any
    sum(c_i * b_i) + tau with full support, 0 < c_i < ord(b_i)/2, and tau a
    word in the basis order-2 generators (stripping an order-2 element is
    orientation-free, and -tau = tau, so both signs of the pin strip to the
    true corners).  Only runs after the standard pin fails on every basis;
    returns ('ok', longs, twos, pin_label, ls_rig) or None.

    Deterministic enumeration order: bases in DFS order; per basis
    coefficient vectors in ascending lexicographic order, tau words in binary
    order (empty first); the all-ones/empty combination is skipped (that is
    the standard pin, already tried)."""
    for basis in _iter_bases_desc(d, cands, order):
        longs = sorted([c for c in basis if c[2] > 2], key=lambda c: c[2], reverse=True)
        twos = [c for c in basis if c[2] == 2]
        if len(longs) < 2:
            continue                        # nothing to pin -- pass 1 handled it
        used = {c[1] for c in basis}
        cs_ranges = [range(1, max((o + 1) // 2, 2)) for _, _, o in longs]
        tau_subsets = list(itertools.chain.from_iterable(
            itertools.combinations(range(len(twos)), k) for k in range(len(twos) + 1)))
        for cs in itertools.product(*cs_ranges):
            for tsub in tau_subsets:
                if all(c == 1 for c in cs) and not tsub:
                    continue                # the standard pin, already tried
                corners = _signed_sum_corners(
                    d, [_pow_desc(dd, c) for (dd, _, _), c in zip(longs, cs)], qf0)
                corners.discard(qf0)
                for ti in tsub:             # shift the corner set by tau
                    tdesc = twos[ti][0]
                    corners = {_desc_qf_pm(q, tdesc)[0] for q in corners}
                for pdesc, pbase, m in cands:
                    if pbase in used:
                        continue
                    if any(qf in corners for qf in _desc_qf_pm(qf0, pdesc) if qf != qf0):
                        pairs = tuple((_pow_desc(dd, c), c)
                                      for (dd, _, _), c in zip(longs, cs))
                        taus = tuple(twos[ti][0] for ti in tsub)
                        pin = ('pin', pdesc, pairs, taus)
                        ls_rig = (tuple(c[0] for c in longs) + (pin,)
                                  + tuple(c[0] for c in twos))
                        return ('ok', longs, twos, pin, ls_rig)
    return None


def disc_rigid_lset_search(d, ls=None):
    """Search for a rigid spanning l-set for the class group of discriminant d.

    ls defaults to modular_prime_pool(): the 15 Atkin primes plus every prime
    with a classical Phi_l in the cache, so the pool grows as new modular
    polynomials are minted.  (Already-cached results are never re-searched by
    the callers, so enlarging the pool only affects fresh searches.)

    Returns a dict.  On success (`success=True`) the key `ls_rig` is a flat
    tuple of primes ready to pass to `ecqf_full_bijection_ord(a, p, ls_rig)`:
        ls_rig = long basis ... , l_sum (if needed), order-2 generators ...
    Other keys: `ls_basis` (independent order>2 gens), `ls_2tors` (independent
    2-torsion gens), `l_sum`, `needs_sum`, `ns` (orders aligned with the
    generators, excluding the pinning prime), and `order` = |A|.

    On failure `success=False`, `message` says which step failed and `best`
    holds the best independent subset found."""
    if ls is None:
        ls = modular_prime_pool()
    qf0 = class_group_id(d)
    A = get_qfs_strict(d)
    order = len(A)
    out = {'d': d, 'order': order, 'success': False, 'message': '',
           'ls_basis': (), 'ls_2tors': (), 'l_sum': None, 'needs_sum': False,
           'ls': (), 'ns': (), 'ls_rig': ()}

    if order == 1:
        out.update(success=True, message='Trivial class group.')
        return out

    cand = _rig_candidates(d, ls, qf0)
    frees = _sib_candidates(d) + _lift_candidates(d, qf0)
    if len(cand) == 0 and len(frees) == 0:
        out['message'] = 'Failed (candidates): no prime in the pool represents a nontrivial element.'
        return out

    # Try bare primes first -- this reproduces the historical search exactly, so
    # discriminants that already succeeded are unchanged.  Fall back to power
    # generators (l, k) when the primes alone admit no independent basis, and to
    # the free generators -- siblings ('sib', q) and lifts ('lift', q) -- when
    # the powers fail too.
    result = None
    if cand:
        result = _rigid_from_candidates(d, qf0, order, _prime_candidates(cand))
        if not (result is not None and result[0] == 'ok'):
            result = _rigid_from_candidates(d, qf0, order, _power_candidates(cand))
    if not (result is not None and result[0] == 'ok') and frees:
        res3 = _rigid_from_candidates(d, qf0, order, _power_candidates(cand) + frees)
        if res3 is not None and (result is None or res3[0] == 'ok'):
            result = res3
    if not (result is not None and result[0] == 'ok'):
        # Last resort: the relaxed pin (any full-support combination of the
        # long generators, coefficients below half the order, plus an optional
        # 2-torsion word) over the richest candidate set.
        res4 = _rigid_relaxed_pin(d, qf0, order, _power_candidates(cand) + frees)
        if res4 is not None:
            result = res4

    if result is None:                              # no independent basis, even with powers
        cand_sorted = sorted(cand, key=lambda l: cand[l], reverse=True)
        best = _best_independent_partial(d, cand_sorted, cand)
        out.update(message='Failed (spanning): the pool has no independent '
                           'generating set for A (even using prime powers and '
                           'the free sibling/lift generators).',
                   best=tuple(best),
                   best_order=len(gen_qfs_d_ls(d, list(best))))
        return out

    if result[0] == 'rigidfail':
        _, longs, twos = result
        out.update(message='Failed (rigid): found a basis with >= 2 generators of '
                           'order > 2, but no candidate represents the sum element '
                           'x* = x_1 + ... + x_m.',
                   best=tuple(c[0] for c in longs) + tuple(c[0] for c in twos),
                   ls_basis=tuple(c[0] for c in longs),
                   ls_2tors=tuple(c[0] for c in twos), needs_sum=True)
        return out

    _, longs, twos, l_sum, ls_rig = result
    out.update(success=True,
               message='Found a rigid spanning set.',
               ls_basis=tuple(c[0] for c in longs),
               ls_2tors=tuple(c[0] for c in twos),
               l_sum=l_sum, needs_sum=(l_sum is not None),
               ls=tuple(c[0] for c in longs) + tuple(c[0] for c in twos),
               ns=tuple(c[2] for c in longs) + tuple(c[2] for c in twos),
               ls_rig=ls_rig)
    return out


def disc_ldata(d, ls=None):
    """disc_rigid_lset_search(d) result, served from qf_ldata.json when available.

    With the default prime pool (ls=None) a precomputed entry is returned if d is
    in qf_ldata, otherwise the search is run live.  Pass an explicit ls to force a
    live search against a custom pool (the cache assumes the default pool)."""
    if ls is None:
        if d in qf_ldata:
            return qf_ldata[d]
        return disc_rigid_lset_search(d)
    return disc_rigid_lset_search(d, ls)


def _velu_cost(a, p, l, walk=False):
    """Velu cost of the prime l for the class of trace a over F_p: the eigenline
    extension degree the data builder will pay, or None if l does not split (no
    horizontal l-isogeny over F_p).  a = 0 is the supersingular case.  With
    walk=False both +-directions are built, so the cost is the LARGER of the two
    degrees; with walk=True (the one-direction cycle walk, velu_nbr_data_ss_walk)
    only the cheap eigenline is used, so the cost is the SMALLER."""
    from nt import frob_ext_degrees
    info = frob_ext_degrees(a, p, l)
    if info['kind'] != 'split':
        return None
    return min(info['degrees']) if walk else max(info['degrees'])

def _velu_rigid_lset_info(a, p, d, deg_caps=(1, 2, 4, 6, 8, 12, 18), prime_cap=4000,
                          walk=False):
    """disc_rigid_lset_search result for discriminant d, drawn from split primes of small
    Velu cost so the eigenline isogenies stay in a manageable extension.  Escalates the cost
    cap until the group is spanned: the default 15-prime pool can force a huge-extension
    generator, so widening to low-cost split primes -- including the cheap ones with small
    eigenvalue order mod l -- avoids it.  The basis generators need a full graph, but the
    sum/pinning prime (info['l_sum']) is only read at the root, so it need not be low-cost.
    A spanning set always exists in principle, so the only obstruction is every generating
    set needing a too-large kernel field (raised as ValueError)."""
    primes = [l for l in primesBetween(2, prime_cap) if l != p]
    cost = {l: _velu_cost(a, p, l, walk=walk) for l in primes}
    for cap in deg_caps:
        pool = tuple(sorted((l for l in primes if cost[l] is not None and cost[l] <= cap),
                            key=lambda l: (cost[l], l)))
        info = disc_rigid_lset_search(d, pool)
        if info.get('success'):
            return info
    raise ValueError(f'no Velu rigid l-set for (a, p) = ({a}, {p}), d = {d} '
                     f'within kernel-field degree {deg_caps[-1]}')

def ss_rigid_lset_info(p, d, **kw):
    """Velu-optimised rigid l-set search for the supersingular class over F_p (trace 0).
    walk=True: the SS data builders walk one eigenline direction per cycle, so a
    generator costs its MIN eigenline degree, not the max."""
    kw.setdefault('walk', True)
    return _velu_rigid_lset_info(0, p, d, **kw)

def ss_rigid_lset(p, d, **kw):
    """The rigid l-set tuple ls_rig for the supersingular class (see ss_rigid_lset_info)."""
    return ss_rigid_lset_info(p, d, **kw)['ls_rig']

def ord_rigid_lset_info(a, p, d=None, **kw):
    """Velu-optimised rigid l-set search for the ordinary class (a, p) -- the fallback used
    when the 15 Atkin primes do not generate the class group of d = a^2 - 4p."""
    if d is None:
        d = a * a - 4 * p
    return _velu_rigid_lset_info(a, p, d, **kw)

def ord_rigid_lset(a, p, d=None, **kw):
    """The rigid l-set tuple ls_rig for the ordinary class (a, p) (see ord_rigid_lset_info)."""
    return ord_rigid_lset_info(a, p, d, **kw)['ls_rig']

def ord_lset(a, p):
    """A rigid l-set for the ordinary class (a, p): the 15 Atkin primes when they span the
    class group (so the curve side can use Atkin modular polynomials), otherwise a
    Velu-optimised set (curve side computed with Velu).  Mirrors the dispatch the pipeline
    uses; `ecqf_full_bijection_ord` falls back to this when no l-set is supplied."""
    d = a * a - 4 * p
    info = disc_ldata(d)                                   # default = the 15 Atkin primes
    if info.get('success'):
        return info['ls_rig']
    return ord_rigid_lset(a, p, d)

def _ss_root_neighbours(p, root, l):
    """The (de-duplicated) horizontal l-isogeny neighbour signatures OF THE ROOT only.
    compute_bijection_zn reads the sum/pinning prime's data only at the root, so this is
    all that is needed for it -- no full graph, hence no cost even for a high-degree l.
    (When the root is the j=1728 curve the two neighbours share a j-invariant.)"""
    from velu import velu_l_isog_codomain
    from ecfp import signature
    f, g = js_to_fg(root, p)
    nbrs = []
    for direction in (0, 1):
        res = velu_l_isog_codomain(f, g, l, p, 0, direction=direction)
        if res['status'] == 'ok' and isinstance(res['f'], int) and isinstance(res['g'], int):
            nbrs.append((res['j'] % p, signature(res['f'], res['g'], p)))
    return list(dict.fromkeys(nbrs))

def _ss_sib_nbrs(p, sig, q):
    """The other members of sig's q-descent fiber on the SS side: the floor
    signatures sharing sig's surface q-parent ("ascend, descend a different
    leg").  Entirely rational 2-isogenies over F_p -- no extension fields --
    so a sib direction is essentially free, mirroring qf_q_siblings."""
    from velu import two_isogeny_sigs
    if q != 2:
        raise NotImplementedError('SS sibling data exists only for q = 2')
    parent = [nb for nb in two_isogeny_sigs(*js_to_fg(sig, p), p)
              if ss_level(nb, p) == 'surface']
    kids = [nb for nb in two_isogeny_sigs(*js_to_fg(parent[0], p), p)
            if ss_level(nb, p) == 'floor']
    return [s for s in dict.fromkeys(kids) if s != sig]

def _ss_sum_root_data(p, root, desc):
    """Root-only neighbour data for a sum/pinning direction: the +-1-step of a
    plain prime, the fiber-mates of a sibling descriptor ('sib', q), or the
    +-k-step of a power (l, k).  When the two eigenline degrees differ, the
    root's whole l-cycle is walked in the cheap direction (velu_walk_cycle) and
    both neighbours are read off it -- n cheap isogenies instead of paying the
    expensive eigenline; with equal degrees the +-steps are computed directly."""
    if _is_sib(desc):
        return _ss_sib_nbrs(p, root, desc[1])
    l, k = desc if isinstance(desc, tuple) else (desc, 1)
    from nt import frob_ext_degrees
    degs = frob_ext_degrees(0, p, l)['degrees']
    if degs[0] != degs[1]:
        from velu import velu_walk_cycle
        cycle, d = velu_walk_cycle(p, l, root)
        n = len(cycle)
        pair = [cycle[k % n], cycle[-k % n]]
        if d == 1:
            pair.reverse()
        return list(dict.fromkeys(pair))
    if k == 1:
        return _ss_root_neighbours(p, root, l)
    out = []
    for first in _ss_root_neighbours(p, root, l):
        prev, cur = root, first
        for _ in range(k - 1):
            step = [s for s in _ss_root_neighbours(p, cur, l) if s != prev] or [prev]
            prev, cur = cur, step[0]
        out.append(cur)
    return list(dict.fromkeys(out))

def _ss_zn_to_sig(p, d, sigs, root):
    """Signature-side Z/n labelling for the SS class on disc d, rooted at `root`, with
    the sum/pinning prime computed only at the root (the only place it is read)."""
    info = ss_rigid_lset_info(p, d)
    ls, l_sum = info['ls_rig'], info.get('l_sum')
    basis = [x for x in ls if x != l_sum]
    sig_nbr = ss_horizontal_nbr_data(p, basis, sigs)
    if l_sum is not None:
        if _is_pin(l_sum):
            # ('pin', pdesc, pairs, taus): pdesc's graph is read only at the
            # root; each pairs key (pk, c) is the c-step of a basis generator,
            # so its graph is the c-step of data already built (no new Velu).
            _, pdesc, pairs, taus = l_sum
            for pk, c in pairs:
                if pk not in sig_nbr:
                    l, m = pk
                    dd = l if m == c else (l, m // c)
                    sig_nbr[pk] = _kstep_nbrdata(sig_nbr[dd], c)
            sig_nbr[l_sum] = {root: _ss_sum_root_data(p, root, pdesc)}
        else:
            sig_nbr[l_sum] = {root: _ss_sum_root_data(p, root, l_sum)}
    return ls, compute_bijection_zn(ls, sig_nbr, root)


        #########################################
        # Step 2: Collect Relevant Isogeny Data #
        #########################################

# Quadratic form/lattice data

def qf_isog_data(d,ls):
    qfs = get_qfs_strict(d)
    # De-duplicate each neighbour list: for an order-2 element x_l the +x_l and
    # -x_l isogenies land on the same class, and qf_isogs_hor returns it twice.
    # Collapsing the duplicate lets l2_split correctly read such an l as a
    # length-1 (order-2) direction rather than mistaking it for a long cycle.
    # A descriptor in ls may also be a power (l, k): its data is the k-step of
    # the base l-graph (no new isogeny data).  A sibling descriptor ('sib', q)
    # gets the q-descent fibers (each class's other q-siblings).
    out = {}
    base = {}
    for desc in ls:
        if _is_pin(desc):
            _, pdesc, pairs, taus = desc
            inner = qf_isog_data(d, (pdesc,) + tuple(pk for pk, c in pairs))
            out[desc] = inner[pdesc]            # the pin label carries pdesc's graph
            for pk, c in pairs:
                out[pk] = inner[pk]
            continue
        if _is_powkey(desc):
            _, bdesc, c = desc
            inner = qf_isog_data(d, (bdesc,))
            out[desc] = _kstep_nbrdata(inner[bdesc], c)
            continue
        if _is_sib(desc):
            out[desc] = {qf: [x for x in qf_q_siblings(qf, desc[1]) if x != qf]
                         for qf in qfs}
            continue
        if _is_lift(desc):
            out[desc] = {qf: qf_lift_nbrs(qf, desc[1]) for qf in qfs}
            continue
        l = _desc_base(desc)
        if l not in base:
            base[l] = {qf: list(dict.fromkeys(qf_isogs_hor(qf, l))) for qf in qfs}
        out[desc] = _kstep_nbrdata(base[l], desc[1]) if isinstance(desc, tuple) else base[l]
    return out


######################
# ECFP Ordinary Data #
######################


def js_to_rabs(js,p):
    ab_to_js = {}
    for i,j1 in enumerate(js):
        for j2 in js[i:]:
            ab_to_js[((-j1-j2)%p,(j1*j2)%p)] = (j1,j2)
    return ab_to_js

@lru_cache(maxsize=32)
def _atk_rab_counts(l, p):
    """{eval_atk(x): multiplicity} over x in F_p.  Class-independent, so the O(p)
    Atkin scan is shared by every isogeny class at p instead of being redone per
    class; the per-class X1 cost drops to O(h^2) dictionary lookups."""
    counts = {}
    for x in range(p):
        evx = eval_atk(x, l, p)
        counts[evx] = counts.get(evx, 0) + 1
    return counts


def ecfp_nbr_data_ord_X1(ap,l,jdata = {}):
    """In-class l-isogeny neighbour data {j: [neighbours, with multiplicity]}.

    Dispatch by modular-polynomial format: the 15 genus-0 primes use the Atkin
    polynomial scan; any other prime with a classical Phi_l in the cache uses
    modpoly_roots_among (roots of Phi_l(j, Y) among the class's j's).  The two
    conventions agree exactly EXCEPT at rows j in {0, 1728} (classical root
    multiplicity carries the extra-automorphism factor 3 resp. 2) and at
    self-loops (the Atkin pair scan appends both directions: 2m vs m).  Both
    consumers are insensitive to the difference: the horizontal path
    (ecfp_nbr_data_ord) de-duplicates rows, and the vertical leaf/ancestor
    logic never treats 0 or 1728 as a leaf or mid-tree vertex -- their CM
    discs are fundamental, so they always sit on the crater."""
    a,p = ap
    js = jdata['js'] if len(jdata) > 0 else trfr_to_js(a,p)
    if l in atkin_polys_dict:
        rabs = jdata['rabs'] if len(jdata) > 0 else js_to_rabs(js,p)
        counts = _atk_rab_counts(l, p)
        nbrdata = {j:[] for j in js}
        for evx, (j1, j2) in rabs.items():
            for _ in range(counts.get(evx, 0)):
                nbrdata[j1].append(j2)
                nbrdata[j2].append(j1)
        return nbrdata
    if l in _modpoly_cache:
        nbrdata = {j:[] for j in js}
        for j1 in nbrdata:
            for j2, m in modpoly_roots_among(j1, l, p, js).items():
                nbrdata[j1] += [j2] * m
        return nbrdata
    raise ValueError(f'no modular polynomial for l = {l}: not one of the 15 '
                     f'Atkin primes and no classical Phi_l in the cache')


def ecfp_nbr_data_ord(ap,ls):
    a,p = ap
    js = trfr_to_js(a,p)
    rabs = js_to_rabs(js,p)
    jspc = {'js':js,'rabs':rabs}
    # De-duplicate neighbour lists (an order-2 direction lists the same j twice),
    # mirroring qf_isog_data, so l2_split reads order-2 vs long the same way on the
    # j side as on the qf side.  The vertical/ancestor code keeps multiplicities by
    # calling ecfp_nbr_data_ord_X1 directly, so it is unaffected by this.  A power
    # descriptor (l, k) is the k-step of the base l-graph -- same data.
    # Dispatch per prime: any prime with a modular polynomial (Atkin format or a
    # classical Phi_l in the cache) goes through the X1 scan; anything else gets
    # its horizontal graph from Velu (velu.velu_nbr_data_ord).  A sibling
    # descriptor ('sib', q) reads the q-descent fibers off the q-ancestor data:
    # js sharing a q-parent are each other's siblings (kernel-coset mates).
    out = {}
    base = {}
    has_modpoly = set(ssprimes) | set(_modpoly_cache)
    for desc in ls:
        if _is_pin(desc):
            _, pdesc, pairs, taus = desc
            inner = ecfp_nbr_data_ord(ap, (pdesc,) + tuple(pk for pk, c in pairs))
            out[desc] = inner[pdesc]            # the pin label carries pdesc's graph
            for pk, c in pairs:
                out[pk] = inner[pk]
            continue
        if _is_powkey(desc):
            _, bdesc, c = desc
            inner = ecfp_nbr_data_ord(ap, (bdesc,))
            out[desc] = _kstep_nbrdata(inner[bdesc], c)
            continue
        if _is_sib(desc) or _is_lift(desc):
            if desc not in base:
                q = desc[1]
                nbrs_q = ecfp_nbr_data_ord_X1(ap, q, jspc)
                anc = tree_edges_to_ancestors(nbrs_q)
                # Fibers over the FLOOR leaves only (the walk's vertex set):
                # tree_edges_to_ancestors can hand crater vertices spurious
                # ancestors when the crater cycle is short, which would pollute
                # the leaf fibers.
                leaves = [j for j in nbrs_q if len(nbrs_q[j]) == 1 and j != 0]
                fibers = {}
                for j in leaves:
                    fibers.setdefault(anc[j], []).append(j)
                if _is_sib(desc):
                    sibmap = {}
                    for grp in fibers.values():
                        for j in grp:
                            sibmap[j] = [j2 for j2 in grp if j2 != j]
                    base[desc] = sibmap
                else:
                    # lifted x_q: j1 -> j2 iff their q-parents are horizontal
                    # q-neighbours; kernel-1, so each parent has one child.
                    leafset = set(leaves)
                    child = {P: grp[0] for P, grp in fibers.items()}
                    liftmap = {}
                    for j in leaves:
                        ups = [x for x in dict.fromkeys(nbrs_q[anc[j]])
                               if x not in leafset]
                        liftmap[j] = list(dict.fromkeys(child[P] for P in ups
                                                        if P in child))
                    base[desc] = liftmap
            out[desc] = base[desc]
            continue
        l = _desc_base(desc)
        if l not in base:
            if l in has_modpoly:
                nbr = ecfp_nbr_data_ord_X1(ap, l, jspc)
                base[l] = {j: list(dict.fromkeys(nbr[j])) for j in nbr}
            else:
                from velu import velu_nbr_data_ord
                base[l] = velu_nbr_data_ord(ap, l, js)
        out[desc] = _kstep_nbrdata(base[l], desc[1]) if isinstance(desc, tuple) else base[l]
    return out

def tree_edges_to_ancestors(nbrdata):
    leaves = [v for v in nbrdata if len(nbrdata[v])==1 and v != 0]
    anc_data = {}
    nextbatch = []
    for v0 in leaves:
        v1 = nbrdata[v0][0]
        anc_data[v0]=v1
        nextbatch.append(v1)
    while len(nextbatch)>0:
        currentbatch = nextbatch
        nextbatch = []
        for v in currentbatch:
            v_ancs = [v1 for v1 in nbrdata[v] if v1 not in anc_data]
            if len(v_ancs)==1:
                anc_data[v]=v_ancs[0]
                nextbatch+=v_ancs
    return anc_data



def get_ancestor_data_ord(ap):
    a,p = ap
    d,c = discfac(a*a-4*p)
    js = trfr_to_js(a,p)
    if c == 1:
        return {'ancestor_data':{},'leaves':js,'js_all':js}
    elif c == 2:
        if d == -4:
            return {'ancestor_data':{287496%p:1728%p},'leaves':[287496%p],'js_all':js}
        leaves = [j for j in js if quad_rec(j-1728,p)==-1]
        ancs = {}
        nbrs = ecfp_nbr_data_ord_X1(ap,2)
        for j in leaves:
            ancs[j]=nbrs[j][0]
        return {'ancestor_data':{2:ancs},'leaves':leaves,'js_all':js}
    ls = [l for l in primesBetween(1,c+1) if c%l ==0]
    anc_data = {}
    leaf_cands = [j for j in js if j!= 0 and (j-1728)%p !=0]
    for l in ls:
        if (c == l and d in heeg_js
                and l not in atkin_polys_dict and l not in _modpoly_cache):
            # Trivial-volcano shortcut: the conductor is exactly one prime l
            # with no modular polynomial, over a class-number-1 fundamental
            # disc d.  The l-volcano then has two levels -- a SINGLE crater
            # vertex, the known CM point j(d) (root of the linear H_d) -- so
            # no Phi_l is needed: everything else is a leaf, and every leaf's
            # l-ancestor is the crater point.  (vert_isog_ext labels the
            # crater via qf_parents, forced to the principal form of d.)
            # Deeper/mixed conductors (e.g. -12 l^2) are NOT covered: their
            # l-crater has several vertices and ascent needs identification.
            crater_j = heeg_js[d] % p
            assert crater_j in set(js), f'CM point j({d}) not in class {ap}'
            anc_data[l] = {j: crater_j for j in js if j != crater_j}
            leaf_cands = [j for j in leaf_cands if j != crater_j]
        else:
            nbrs_l = ecfp_nbr_data_ord_X1((a,p),l)
            anc_data[l]=tree_edges_to_ancestors(nbrs_l)
            leaf_cands = [j for j in leaf_cands if len(nbrs_l[j])==1]
    return {'ancestor_data':anc_data,'leaves':leaf_cands,'js_all':js}

def zn_ecqf_bij(a:int,p:int,ls:tuple[int]):
    anc_data = get_ancestor_data_ord((a,p))
    j0 = anc_data['leaves'][0]
    isog_data_horz_fp = ecfp_nbr_data_ord((a,p),ls)
    zn_to_j = compute_bijection_zn(ls,isog_data_horz_fp,j0)
    zn_to_qf = compute_bijection_zn(ls,qf_isog_data(a*a-4*p,ls),class_group_id(a*a-4*p))
    assert len(zn_to_j)==len(zn_to_qf)
    assert len({t for t in zn_to_qf if t not in zn_to_j})== 0
    return {t:(zn_to_j[t],zn_to_qf[t]) for t in zn_to_qf}


def vert_isog_ext(j_to_qf:dict,vertical_iso_data:dict)->dict:
    for l in vertical_iso_data['ancestor_data']:
        ancsl = vertical_iso_data['ancestor_data'][l]
        nextbatch = [j for j in j_to_qf if j in ancsl and ancsl[j] not in j_to_qf]
        while len(nextbatch)>0:
            currentbatch = nextbatch.copy()
            nextbatch = []
            for j0 in currentbatch:
                j1 = ancsl[j0]
                if j1 not in j_to_qf:
                    qf0 = j_to_qf[j0]
                    qf1 = qf_parents(qf0,l)[0]
                    j_to_qf[j1] = qf1
                    if j1 in ancsl:
                        nextbatch.append(j1)
    return j_to_qf

        #########################################
        # Orientation canonicalization          #
        #########################################
#
# A rigid l-set pins the cube labelling up to the single global automorphism
# x -> -x (complex conjugation / class-group inversion); see the rigid l-set
# notes above.  So each of the two labellings (qf side and j side) has exactly
# one binary orientation freedom.  We fix both with deterministic, intrinsic
# rules so the end-to-end bijection is reproducible and -- crucially -- the
# qf-side labelling is a function of d alone, hence consistent across every
# (a,p) in the same discriminant (and so across isogeny classes).
#
# Convention:
#   * qf side: the first long generator e_1 = (1,0,...,0) maps to the reduced
#     form with b > 0.  (e_1's image has order > 2, so it is non-ambiguous and
#     b != 0; its inverse is the b < 0 form.)
#   * j side: the root's +1 neighbour in coordinate 0 is the numerically
#     smaller j-invariant.
# When there is no generator of order > 2 the group is 2-torsion, every form is
# ambiguous, the j-invariants are all "real", and conjugation is trivial -- both
# rules then leave the labelling untouched.

def _labelling_orders(bij):
    keys = list(bij)
    n = len(keys[0])
    return tuple(1 + max(k[i] for k in keys) for i in range(n))

def _neg_tuple(t,orders):
    return tuple((-ti) % n for ti,n in zip(t,orders))

def canonicalize_qf_labelling(zn_to_qf):
    """Pin the global x->-x freedom on the qf side (b>0 at the first long gen)."""
    if not zn_to_qf:
        return zn_to_qf
    n = len(next(iter(zn_to_qf)))
    if n == 0:
        return zn_to_qf
    e1 = (1,)+(0,)*(n-1)
    q1 = zn_to_qf.get(e1)
    if q1 is None or class_group_inv(q1) == q1:   # ambiguous => no order>2 coord
        return zn_to_qf
    if q1[1] < 0:
        return {t:class_group_inv(v) for t,v in zn_to_qf.items()}
    return zn_to_qf

def canonicalize_j_labelling(zn_to_j):
    """Pin the same global freedom on the j side (smaller j at root's +1 nbr)."""
    if not zn_to_j:
        return zn_to_j
    orders = _labelling_orders(zn_to_j)
    if len(orders) == 0 or orders[0] <= 2:        # coord 0 is 2-torsion => trivial
        return zn_to_j
    n = len(orders)
    e1  = (1,)+(0,)*(n-1)
    em1 = (orders[0]-1,)+(0,)*(n-1)
    if zn_to_j[e1] > zn_to_j[em1]:
        return {_neg_tuple(t,orders):v for t,v in zn_to_j.items()}
    return zn_to_j


def ecqf_full_bijection_ord(a:int,p:int,ls:tuple[int]=None,zn_to_qf=None):
    # zn_to_qf is the (a,p)-independent labelling tuple -> qf of the class group.
    # It depends only on (d, ls), so it can be precomputed once per discriminant
    # and reused for every (a,p); pass it in to skip the expensive recompute.
    # ls defaults to ord_lset(a, p): the 15 Atkin primes when they span, else Velu.
    if a == 0:
        return ecqf_full_bijection_ss(p)
    d = a*a-4*p
    hd = small_bij_check(d)
    if len(hd)>0:
        return {j%p:hd[j] for j in hd}
    if ls is None:
        ls = ord_lset(a, p)
    vert_iso_data_fp = get_ancestor_data_ord((a,p))
    j0 = vert_iso_data_fp['leaves'][0]
    isog_data_horz_fp = ecfp_nbr_data_ord((a,p),ls)
    zn_to_j = canonicalize_j_labelling(compute_bijection_zn(ls,isog_data_horz_fp,j0))
    if zn_to_qf is None:
        zn_to_qf = compute_bijection_zn(ls,qf_isog_data(d,ls),class_group_id(d))
    zn_to_qf = canonicalize_qf_labelling(zn_to_qf)   # idempotent if already canonical
    assert len(zn_to_j)==len(zn_to_qf)
    assert len({t for t in zn_to_qf if t not in zn_to_j})== 0
    return vert_isog_ext(j_to_qf=compdiv_dics(zn_to_j,zn_to_qf),vertical_iso_data=vert_iso_data_fp)


#########################
# Supersingular over F_p #
#########################
# Objects are SIGNATURES (j, s), s = +-1 (an F_p-iso class = a curve or its twist);
# the isogeny class is an oriented CM torsor under the class group of an order in
# Q(sqrt(-p)).  Trace is 0, so neighbour data comes from Velu (velu_nbr_data_ss),
# not the Atkin polynomials.  Case A (p = 1 mod 4) has conductor 1 -- every curve is
# maximal, Z[sqrt(-p)], disc -4p -- so there is no volcano and the assembler is a
# direct signature analogue of zn_ecqf_bij with no ancestor/vertical machinery.

def ss_signatures(p):
    """Signature set of the whole supersingular isogeny class over F_p: both signs
    of every F_p-rational SS j-invariant (curve and quadratic twist).  Spans both
    orders of Q(sqrt(-p)) when p = 3 mod 4."""
    return [(j, s) for j in trfr_to_js(0, p) for s in (1, -1)]

def ss_level(sig, p):
    """'surface' (maximal order O_K, full 2-torsion) or 'floor' (Z[sqrt(-p)], single
    order-2 point) for a signature, via the j-1728 criterion (p = 3 mod 4 only)."""
    j, s = sig
    if (j - 1728) % p == 0:
        return 'surface' if s == -1 else 'floor'      # (1728,-1) full 2-tors; (1728,+1) single
    return 'surface' if quad_rec((j - 1728) % p, p) == 1 else 'floor'

def ss_horizontal_nbr_data(p, ls, sigs):
    """Horizontal isogeny neighbour data {desc: {sig:[nbr sigs]}} over a SET of
    same-level signatures.  Odd l via the eigenline Velu (velu_nbr_data_ss); l = 2 via
    the rational 2-torsion 2-isogenies (two_isogeny_sigs), kept to same-level codomains
    (the horizontal 2-isogenies, which exist on the surface when 2 splits in O_K)."""
    from velu import velu_nbr_data_ss_walk, two_isogeny_sigs
    sset = set(sigs)
    out, base = {}, {}
    for desc in ls:
        if _is_pin(desc):
            # No root known here, so pdesc gets its full graph (the root-only
            # shortcut lives in _ss_zn_to_sig); mirrors qf_isog_data.
            _, pdesc, pairs, taus = desc
            inner = ss_horizontal_nbr_data(p, (pdesc,) + tuple(pk for pk, c in pairs), sigs)
            out[desc] = inner[pdesc]            # the pin label carries pdesc's graph
            for pk, c in pairs:
                out[pk] = inner[pk]
            continue
        if _is_sib(desc):
            out[desc] = {sig: _ss_sib_nbrs(p, sig, desc[1]) for sig in sigs}
            continue
        if _is_lift(desc) or _is_powkey(desc):
            raise NotImplementedError(
                f'SS signature side has no data builder for descriptor {desc!r}')
        l = _desc_base(desc)
        if l not in base:
            if l == 2:
                base[2] = {sig: list(dict.fromkeys(
                    nb for nb in two_isogeny_sigs(*js_to_fg(sig, p), p) if nb in sset))
                    for sig in sigs}
            else:
                base[l] = velu_nbr_data_ss_walk(p, l, sigs)
        out[desc] = _kstep_nbrdata(base[l], desc[1]) if isinstance(desc, tuple) else base[l]
    return out

def canonicalize_sig_labelling(zn_to_sig):
    """Pin the global x->-x freedom on the signature side, as canonicalize_j_labelling
    does on the j side: the root's +1 neighbour in coord 0 is the lexicographically
    smaller signature (j, then s)."""
    if not zn_to_sig:
        return zn_to_sig
    orders = _labelling_orders(zn_to_sig)
    if len(orders) == 0 or orders[0] <= 2:
        return zn_to_sig
    n = len(orders)
    e1  = (1,) + (0,) * (n - 1)
    em1 = (orders[0] - 1,) + (0,) * (n - 1)
    if zn_to_sig[e1] > zn_to_sig[em1]:
        return {_neg_tuple(t, orders): v for t, v in zn_to_sig.items()}
    return zn_to_sig

def ecqf_full_bijection_ss(p, ls=None):
    """Supersingular signature <-> qf bijection over F_p, from scratch via Velu.
    Case A (p = 1 mod 4): a signature analogue of zn_ecqf_bij (no volcano).
    Case B (p = 7 mod 8): surface bijection + descend to the floor by 2-isogeny.
    Case C (p = 3 mod 8): floor bijection + ascend to the surface.
    `ls` overrides the (SS-cost-aware) rigid l-set; otherwise ss_rigid_lset picks one."""
    if p % 4 == 1:
        d = -4 * p
        sigs = ss_signatures(p)
        root = min(sigs)
        if ls is None:
            ls, zn_to_sig = _ss_zn_to_sig(p, d, sigs, root)
        else:
            zn_to_sig = compute_bijection_zn(ls, ss_horizontal_nbr_data(p, ls, sigs), root)
        zn_to_sig = canonicalize_sig_labelling(zn_to_sig)
        zn_to_qf = canonicalize_qf_labelling(
            compute_bijection_zn(ls, qf_isog_data(d, ls), class_group_id(d)))
        assert len(zn_to_sig) == len(zn_to_qf)
        assert len({t for t in zn_to_qf if t not in zn_to_sig}) == 0
        return compdiv_dics(zn_to_sig, zn_to_qf)
    if p % 8 == 7:
        return ecqf_full_bijection_ss_B(p, ls)
    return ecqf_full_bijection_ss_C(p, ls)


def ecqf_full_bijection_ss_C(p, ls_floor=None):
    """Case C (p = 3 mod 8): 2 is inert in O_K, so h(-4p) = 3*h(-p) and the surface has
    no horizontal 2-isogenies.  The FLOOR is the order Z[sqrt(-p)] (disc -4p) carrying
    the full class group, so its bijection is a case-A-style assembly on disc -4p (the
    rigid l-set search already handles a non-cyclic Cl(-4p) via its sum/pinning prime),
    rooted at (1728, +1).  Then ascend each floor node to its unique surface parent."""
    from velu import two_isogeny_sigs
    d_floor = -4 * p
    floor = [s for s in ss_signatures(p) if ss_level(s, p) == 'floor']
    root = (1728 % p, 1)
    if ls_floor is None:
        ls_floor, zn_to_sig = _ss_zn_to_sig(p, d_floor, floor, root)
    else:
        zn_to_sig = compute_bijection_zn(ls_floor, ss_horizontal_nbr_data(p, ls_floor, floor), root)
    zn_to_qf = compute_bijection_zn(ls_floor, qf_isog_data(d_floor, ls_floor),
                                    class_group_id(d_floor))
    floor_to_qf = compdiv_dics(zn_to_sig, zn_to_qf)
    # Ascend each floor node to its unique surface parent (3 floor nodes share a parent;
    # the assignments agree, so dedup is automatic).
    sig_to_qf = dict(floor_to_qf)
    for fs, fq in floor_to_qf.items():
        par_sig = [nb for nb in two_isogeny_sigs(*js_to_fg(fs, p), p) if ss_level(nb, p) == 'surface']
        par_qf = [qf_mod_gamma(q) for q in qf_isogs(fq, 2) if qf_disc(q) == -p]
        sig_to_qf[par_sig[0]] = par_qf[0]
    return sig_to_qf


def ecqf_full_bijection_ss_B(p, ls_surf=None):
    """Case B (p = 7 mod 8): h(-p) = h(-4p), depth-1 volcano, each surface node has
    exactly one descending 2-isogeny.  Build the bijection on the SURFACE (maximal
    order O_K, disc -p) -- 2 splits there so the horizontal 2-isogeny cycle can carry
    it -- rooted at (1728, -1), then extend to the floor (disc -4p) by mapping each
    surface node's unique descending 2-isogeny to the corresponding qf descent."""
    from velu import two_isogeny_sigs
    d_surf, d_floor = -p, -4 * p
    sigs = ss_signatures(p)
    surf = [s for s in sigs if ss_level(s, p) == 'surface']
    root = (1728 % p, -1)
    if ls_surf is None:
        ls_surf, zn_to_sig = _ss_zn_to_sig(p, d_surf, surf, root)
    else:
        zn_to_sig = compute_bijection_zn(ls_surf, ss_horizontal_nbr_data(p, ls_surf, surf), root)
    zn_to_qf = canonicalize_qf_labelling(
        compute_bijection_zn(ls_surf, qf_isog_data(d_surf, ls_surf), class_group_id(d_surf)))
    surf_to_qf = compdiv_dics(zn_to_sig, zn_to_qf)
    # Descend: each surface node's unique floor child <-> the qf form one order down.
    sig_to_qf = dict(surf_to_qf)
    for ssig, sqf in surf_to_qf.items():
        child = [nb for nb in two_isogeny_sigs(*js_to_fg(ssig, p), p)
                 if ss_level(nb, p) == 'floor']
        below = [qf_mod_gamma(q) for q in qf_isogs(sqf, 2) if qf_disc(q) == d_floor]
        sig_to_qf[child[0]] = below[0]
    return sig_to_qf

