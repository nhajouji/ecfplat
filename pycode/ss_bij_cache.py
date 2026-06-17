"""
Populate / update data/ecqf_ss_pcbij_velu_4_1024.json -- the supersingular
signature<->lattice bijection recomputed FROM SCRATCH with the Velu pipeline
(ecqf_bij.ecqf_full_bijection_ss), kept SEPARATE from the original Sage-computed
files (ecqf_ss_pcbij_4_1024_INC.json / ssfp_pc_bij.json).

Format matches the Sage files: {str(p): {"(j, s)": [a, b, c]}}, j-invariant j in
F_p and signature s = +-1, value a reduced binary quadratic form.

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
import argparse

from nt import primesBetween
from ecqf_bij import ecqf_full_bijection_ss, _DATA_DIR
from qfs import qf_mod_gamma

DEFAULT_PATH = _DATA_DIR / 'ecqf_ss_pcbij_velu_4_1024.json'
LISTFORM_PATH = _DATA_DIR / 'ssfp_pc_bij_velu.json'


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
    """Derive the list-form file [[[j, s], [a, b, c]], ...] (mirroring ssfp_pc_bij.json)
    from the dict-form bijection file, so the two stay consistent."""
    data = load(src)
    out = {p: [[list(ast.literal_eval(k)), v] for k, v in entry.items()]
           for p, entry in data.items()}
    save(out, dst)
    return out


def populate(pmin=5, pmax=1024, path=DEFAULT_PATH, force=False, save_every=20, verbose=True):
    data = load(path)
    primes = [p for p in primesBetween(max(pmin, 5), pmax)]
    done = skipped = failed = 0
    for p in primes:
        key = str(p)
        if (not force) and key in data:
            skipped += 1
            continue
        try:
            t0 = time.time()
            data[key] = bijection_entry(p)
            done += 1
            if verbose:
                print(f'p={p:>4}: {len(data[key])} signatures  ({time.time() - t0:.1f}s)', flush=True)
        except Exception as e:                     # keep going; report at the end
            failed += 1
            print(f'p={p:>4}: FAILED {type(e).__name__}: {e}', flush=True)
        if done and done % save_every == 0:
            save(data, path)
    save(data, path)
    write_list_form(path)                          # keep the list-form file in sync
    if verbose:
        print(f'done={done} skipped={skipped} failed={failed} total={len(data)} -> {path}')
    return data


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--min', type=int, default=5, dest='pmin')
    ap.add_argument('--max', type=int, default=1024, dest='pmax')
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    populate(pmin=args.pmin, pmax=args.pmax, force=args.force)
