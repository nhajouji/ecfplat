"""
Populate / update data/qf_ldata.json -- precomputed disc_rigid_lset_search output.

This is the lightweight l-set analog of rigid_cache.py (which manages the heavier
per-discriminant bijection cache, rigid_lset_cache.json).  qf_ldata.json is a flat
JSON map {str(d): disc_rigid_lset_search(d)} that ecqf_bij.py loads at import into the
module dict `qf_ldata` (read it via ecqf_bij.disc_ldata).

CLI
---
    python ldata_cache.py --min -4096            # fill -3 .. -4096, skip what's present
    python ldata_cache.py --min -8192            # widen the range; only fills new d
    python ldata_cache.py --min -4096 --force    # recompute every entry
Re-running is incremental (already-present discriminants are skipped unless --force),
so this doubles as the "update later / more compute" script.  After regenerating,
re-import ecqf_bij to pick up the new data.
"""

import os
import json

from ecqf_bij import disc_rigid_lset_search, ssprimes, _DATA_DIR

DEFAULT_PATH = _DATA_DIR / 'qf_ldata.json'


def load_ldata(path=DEFAULT_PATH):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_ldata(data, path=DEFAULT_PATH):
    tmp = str(path) + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(data, f)
    os.replace(tmp, path)                 # atomic on POSIX


def discriminants_in_range(dmin, dmax=-3):
    """Negative discriminants (== 0,1 mod 4) from dmax down to dmin, inclusive."""
    return [d for d in range(dmax, dmin - 1, -1) if d % 4 in (0, 1)]


def populate(dmin, dmax=-3, path=DEFAULT_PATH, pool=ssprimes,
             force=False, save_every=400, verbose=True):
    """Fill qf_ldata.json with disc_rigid_lset_search(d) for every d in [dmin, dmax].

    Skips discriminants already present (unless force=True), so re-running with a
    wider range or after a speed-up only computes the missing ones."""
    data = load_ldata(path)
    todo = discriminants_in_range(dmin, dmax)
    done = skipped = failed = 0
    for d in todo:
        key = str(d)
        if (not force) and key in data:
            skipped += 1
            continue
        try:
            res = disc_rigid_lset_search(d, pool)
            data[key] = res
            if not res.get('success'):
                failed += 1
        except Exception as e:            # never let one d abort the whole run
            data[key] = {'d': d, 'success': False,
                         'message': f'EXC: {type(e).__name__}: {e}'}
            failed += 1
        done += 1
        if save_every and done % save_every == 0:
            save_ldata(data, path)
            if verbose:
                print(f'  ... {done} computed (at d={d}), {skipped} skipped', flush=True)
    save_ldata(data, path)
    if verbose:
        print(f'done: {done} computed, {skipped} already present, {failed} with no '
              f'rigid l-set; file holds {len(data)} discriminants -> {path}', flush=True)
    return data


def _main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        description='Populate / update data/qf_ldata.json (precomputed disc_rigid_lset_search).')
    ap.add_argument('--min', type=int, default=-4096,
                    help='most negative discriminant, inclusive (default -4096)')
    ap.add_argument('--max', type=int, default=-3,
                    help='least negative discriminant, inclusive (default -3)')
    ap.add_argument('--path', default=DEFAULT_PATH, help='qf_ldata JSON path')
    ap.add_argument('--force', action='store_true',
                    help='recompute entries even if already present')
    ap.add_argument('--save-every', type=int, default=400,
                    help='checkpoint the JSON every N computed entries')
    args = ap.parse_args(argv)
    populate(args.min, args.max, path=args.path, force=args.force,
             save_every=args.save_every)


if __name__ == '__main__':
    _main()
