"""
Populate / update data/ecqf_ss_pcbij_velu_4_1024.json -- the supersingular
signature<->lattice bijection computed FROM SCRATCH with the Velu pipeline
(ecqf_bij.ecqf_full_bijection_ss).  This is the live supersingular data the rest of
the code loads (ecqf_tools.ecqf_ss_1K_pc, ecfp.ss_precomputed_dictionary); it
replaced the project's original Sage-computed tables.

Format: {str(p): {"(j, s)": [a, b, c]}}, j-invariant j in F_p and signature
s = +-1, value a reduced binary quadratic form.

CLI
---
    python ss_bij_cache.py                     # fill primes 5..1024, skip present
    python ss_bij_cache.py --max 2048          # widen the range; only new primes
    python ss_bij_cache.py --force             # recompute every prime
    python ss_bij_cache.py --min 500 --max 800 # a sub-range
Re-running is incremental (primes already present are skipped unless --force), and
the file is saved atomically every `save_every` primes so a long run is resumable.
"""

import os
import ast
import json
import time
import random
import argparse

from nt import primesBetween
from ecqf_bij import ecqf_full_bijection_ss, _DATA_DIR
from qfs import qf_mod_gamma

DEFAULT_PATH = _DATA_DIR / 'ecqf_ss_pcbij_velu_4_1024.json'
LISTFORM_PATH = _DATA_DIR / 'ssfp_pc_bij_velu.json'


        ######################
        # Validation battery #
        ######################

# The primes past 1024 have no Sage ground truth (the original notebook is
# lost), so a fresh entry is accepted on internal evidence.  The battery is
# chosen so each leg is independent of the machinery that BUILT the entry
# (Velu eigenline isogenies + the qf rigid-set assembly):
#   * sigs/forms:   the two sides are exactly the expected sets -- signatures
#                   from the trace-0 model count, forms from the reduced-form
#                   enumeration -- and the map is a bijection.
#   * genuine_ss:   every curve really is supersingular: (p+1) * P = O for
#                   random points P on each signature's canonical model (pure
#                   F_p group arithmetic, no isogenies).
#   * phi_edges:    the graph structure is real: at a split prime lv UNUSED by
#                   the construction, every horizontal lv-edge predicted on the
#                   qf side satisfies Phi_lv(j, j') = 0 mod p.  The classical
#                   Phi_lv comes from the CRT/q-expansion pipeline -- a fully
#                   independent code path from Velu.
#   * root:         (p = 3 mod 4 only) the j = 1728 curve sits on the class
#                   group identity, per the root convention.

def _phi_eval(M, x, y, p):
    """Phi_l(x, y) mod p from the dense coefficient matrix M."""
    xp = [1] * len(M)
    for i in range(1, len(M)):
        xp[i] = xp[i - 1] * x % p
    tot = 0
    for i, row in enumerate(M):
        s = 0
        for jj, co in enumerate(row):
            if co:
                s += co * xp[jj]
        tot += s % p * pow(y, i, p)
    return tot % p


def _entry_sig_to_qf(entry):
    """{(j, s): (a, b, c)} from a stored {\"(j, s)\": [a, b, c]} entry."""
    return {ast.literal_eval(k): tuple(v) for k, v in entry.items()}


def validate_entry(p, entry=None, edge_prime=None, n_points=2, seed=0):
    """Validation battery for one prime's bijection entry (stored form or fresh
    from bijection_entry).  Returns {'p', 'ok', 'checks': {name: bool}, ...};
    never raises on a failed check, only on missing data."""
    from ecqf_bij import (ss_signatures, ss_rigid_lset_info, _desc_base,
                          _is_pin, _is_sib, _is_lift)
    from qfs import get_qfs_strict, qf_isogs_hor, class_group_id
    from modularpolynomials import _modpoly_cache
    from ecfp import js_to_fg
    from alg_classes import GF_p
    from velu import embed_fp, random_point, ec_mul
    if entry is None:
        entry = bijection_entry(p)
    sig_to_qf = _entry_sig_to_qf(entry)
    checks, info = {}, {}

    # -- sigs / forms are exactly the expected sets, bijectively
    sigs = set(ss_signatures(p))
    checks['sigs'] = set(sig_to_qf) == sigs
    discs = (-4 * p,) if p % 4 == 1 else (-p, -4 * p)
    expected = {qf_mod_gamma(q) for d in discs for q in get_qfs_strict(d)}
    forms = [qf_mod_gamma(q) for q in sig_to_qf.values()]
    checks['forms'] = set(forms) == expected and len(set(forms)) == len(forms)

    # -- every curve is genuinely supersingular: (p+1) P = O
    rng = random.Random(seed)
    K = GF_p(p)
    ok = True
    for (j, s) in sig_to_qf:
        f, g = js_to_fg((j, s), p)
        a, b = embed_fp(f, K), embed_fp(g, K)
        for _ in range(n_points):
            if ec_mul(p + 1, random_point(a, b, rng), a) is not None:
                ok = False
    checks['genuine_ss'] = ok

    # -- Phi edge check at an unused split classical-modpoly prime
    used = set()
    for d in discs:
        for x in ss_rigid_lset_info(p, d)['ls_rig']:
            if _is_pin(x):
                pd = x[1]
                if not (_is_sib(pd) or _is_lift(pd)):
                    used.add(_desc_base(pd))
                used.update(_desc_base(pk) for pk, c in x[2])
            elif not (_is_sib(x) or _is_lift(x)):
                used.add(_desc_base(x))
    if edge_prime is None:
        cands = [l for l in sorted(_modpoly_cache) if l % 2 and l != p
                 and l not in used and pow(-p % l, (l - 1) // 2, l) == 1]
        edge_prime = cands[0] if cands else None
    info['edge_prime'] = edge_prime
    if edge_prime is None:
        checks['phi_edges'] = None                 # no unused split Phi_l available
    else:
        M = _modpoly_cache[edge_prime]
        qf_to_sig = {q: sig for sig, q in sig_to_qf.items()}
        ok, n_edges = True, 0
        for sig, q in sig_to_qf.items():
            for nb in dict.fromkeys(qf_isogs_hor(q, edge_prime)):
                nb_sig = qf_to_sig.get(qf_mod_gamma(nb))
                n_edges += 1
                if nb_sig is None or _phi_eval(M, sig[0], nb_sig[0], p) != 0:
                    ok = False
        checks['phi_edges'] = ok
        info['n_edges'] = n_edges

    # -- root convention (p = 3 mod 4): j = 1728 on the class group identity
    if p % 4 == 3:
        d0 = -p if p % 8 == 7 else -4 * p
        s0 = -1 if p % 8 == 7 else 1
        root_qf = sig_to_qf.get((1728 % p, s0))
        checks['root'] = (root_qf is not None
                          and qf_mod_gamma(root_qf) == qf_mod_gamma(class_group_id(d0)))

    return {'p': p, 'ok': all(v for v in checks.values() if v is not None),
            'checks': checks, **info}


def load(path=DEFAULT_PATH):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save(data, path=DEFAULT_PATH):
    tmp = str(path) + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(data, f)
    os.replace(tmp, path)                          # atomic on POSIX


def bijection_entry(p):
    """{"(j, s)": [a, b, c]} for the supersingular class over F_p, via Velu."""
    bij = ecqf_full_bijection_ss(p)
    return {str(sig): list(qf_mod_gamma(qf)) for sig, qf in bij.items()}


def write_list_form(src=DEFAULT_PATH, dst=LISTFORM_PATH):
    """Derive the list-form file [[[j, s], [a, b, c]], ...] from the dict-form bijection
    file, so the two stay consistent."""
    data = load(src)
    out = {p: [[list(ast.literal_eval(k)), v] for k, v in entry.items()]
           for p, entry in data.items()}
    save(out, dst)
    return out


def populate(pmin=5, pmax=1024, path=DEFAULT_PATH, force=False, save_every=20,
             verbose=True, validate=False):
    """Fill the bijection cache over [pmin, pmax], skipping present primes.
    With validate=True each fresh entry must pass the validation battery before
    it is stored; failing primes are reported and left out (re-runnable)."""
    data = load(path)
    primes = [p for p in primesBetween(max(pmin, 5), pmax)]
    done = skipped = failed = 0
    failures = []
    for p in primes:
        key = str(p)
        if (not force) and key in data:
            skipped += 1
            continue
        try:
            t0 = time.time()
            entry = bijection_entry(p)
            tag = ''
            if validate:
                res = validate_entry(p, entry)
                if not res['ok']:
                    bad = [k for k, v in res['checks'].items() if v is False]
                    raise AssertionError(f'validation failed: {bad}')
                tag = f'  [valid, edges at l={res["edge_prime"]}]'
            data[key] = entry
            done += 1
            if verbose:
                print(f'p={p:>4}: {len(entry)} signatures  '
                      f'({time.time() - t0:.1f}s){tag}', flush=True)
        except Exception as e:                     # keep going; report at the end
            failed += 1
            failures.append((p, f'{type(e).__name__}: {e}'))
            print(f'p={p:>4}: FAILED {type(e).__name__}: {e}', flush=True)
        if done and done % save_every == 0:
            save(data, path)
    save(data, path)
    write_list_form(path)                          # keep the list-form file in sync
    if verbose:
        print(f'done={done} skipped={skipped} failed={failed} total={len(data)} -> {path}')
        for p, msg in failures:
            print(f'  failed p={p}: {msg}')
    return data


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--min', type=int, default=5, dest='pmin')
    ap.add_argument('--max', type=int, default=1024, dest='pmax')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--validate', action='store_true',
                    help='run the validation battery on each fresh entry')
    args = ap.parse_args()
    populate(pmin=args.pmin, pmax=args.pmax, force=args.force, validate=args.validate)
