"""Closed-form predictions for the edge census (A, B, C, D, Q) and their
verification against the kernel-level model.

    D = all supersingular points of X_0(l) over F_p-bar
    A = F_p-rational ones (Frobenius-fixed)
    B = Hermitian ones (Frobenius = dual)
    C = self-dual ones (dual-fixed)
    Q = (D + C - A - B)/4 = number of quadratic points of X_0(l)^+ on the ss locus
      = g(X_0(pl)^*) - g(X_0(l)^+)   [the grand formula]

Formulas (p, l >= 5 distinct primes; conventions as in the proof notes):
    A = N_rat(p) * (1 + (-4p/l))
    B = (h(-4pl) + [pl = 3 mod 4] h(-pl)) / 2
    C = (h(-4l) + [l = 3 mod 4] h(-l)) * [p nonsplit in Q(sqrt(-l))]
    D = (S(p) - d0 - d1728)(l+1) + d0 o_3(l) + d1728 o_2(l)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from nt import quad_rec, primesBetween
from identities import clgr_size_gen as _h

def h(d):
    """Class number of the order of discriminant d (< 0, = 0 or 1 mod 4)."""
    assert d < 0 and d % 4 in (0, 1), d
    return _h(d)

def S(p):
    return p // 12 + {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]

def N_rat(p):
    """Number of supersingular j in F_p (Deuring / Nakaya Thm 2)."""
    if p % 4 == 1:
        return h(-4 * p) // 2
    hp = h(-p)
    return 2 * hp if p % 8 == 3 else hp

def A_formula(p, l):
    return N_rat(p) * (1 + quad_rec(-4 * p, l))

def B_formula(p, l):
    tot = h(-4 * p * l) + (h(-p * l) if (p * l) % 4 == 3 else 0)
    assert tot % 2 == 0
    return tot // 2

def C_formula(p, l):
    f = -l if l % 4 == 3 else -4 * l
    if quad_rec(f, p) == 1:
        return 0
    return h(-4 * l) + (h(-l) if l % 4 == 3 else 0)

def D_formula(p, l):
    d0 = 1 if p % 3 == 2 else 0
    d1728 = 1 if p % 4 == 3 else 0
    D = (S(p) - d0 - d1728) * (l + 1)
    if d0:
        f3 = 1 + quad_rec(-3, l); D += f3 + (l + 1 - f3) // 3
    if d1728:
        f2 = 1 + quad_rec(-4, l); D += f2 + (l + 1 - f2) // 2
    return D

def census_formula(p, l):
    A, B, C, D = A_formula(p, l), B_formula(p, l), C_formula(p, l), D_formula(p, l)
    num = D + C - A - B
    assert num % 4 == 0, (p, l, A, B, C, D)
    return {'A': A, 'B': B, 'C': C, 'D': D, 'Q': num // 4,
            'points': (D + C) // 2, 'rational': (A + B) // 2}

def genus_prediction(p, l):
    """Q via the grand formula, for cross-checking the census formula."""
    from pell_duality import genus_x0_star, genus_x0_plus
    return genus_x0_star(p, l) - genus_x0_plus(l)

def prediction_table(ls=(5, 7, 11, 19), pmax=128, verify=False, verify_max_l=None,
                     verbose=True):
    """Predicted (A, B, C, D, Q) for all primes p < pmax and l in ls; optional
    kernel-model verification (slow-ish for large l)."""
    rows = []
    for l in ls:
        for p in primesBetween(5, pmax):
            if p == l:
                continue
            f = census_formula(p, l)
            g = genus_prediction(p, l)
            row = {'p': p, 'l': l, **f, 'Q_genus': g, 'consistent': f['Q'] == g}
            if verify and (verify_max_l is None or l <= verify_max_l):
                from pell_duality import involution_profile
                r = involution_profile(p, l)
                row['verified'] = all(r[k] == f[k] for k in 'ABCD') and r['Q_full'] == f['Q']
                row['model'] = {k: r[k] for k in 'ABCD'} | {'Q': r['Q_full']}
            rows.append(row)
    if verbose:
        for l in ls:
            print(f'\n===== l = {l} =====')
            print(f"{'p':>4} {'A':>3} {'B':>3} {'C':>3} {'D':>4} {'Q':>3}  {'Q(genus)':>8}  {'notes'}")
            for r in rows:
                if r['l'] != l:
                    continue
                notes = []
                if not r['consistent']:
                    notes.append('!! census != genus')
                if 'verified' in r:
                    notes.append('model OK' if r['verified'] else f"!! model {r['model']}")
                if p_is_monster(r['p']):
                    notes.append('(Monster p)')
                print(f"{r['p']:>4} {r['A']:>3} {r['B']:>3} {r['C']:>3} {r['D']:>4} {r['Q']:>3}  {r['Q_genus']:>8}  {' '.join(notes)}")
    return rows

MONSTER = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71}
def p_is_monster(p):
    return p in MONSTER

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--verify-max-l', type=int, default=None)
    ap.add_argument('--pmax', type=int, default=128)
    args = ap.parse_args()
    rows = prediction_table(verify=args.verify, verify_max_l=args.verify_max_l, pmax=args.pmax)
    bad = [r for r in rows if not r['consistent'] or r.get('verified') is False]
    print(f'\n{len(rows)} rows; census/genus inconsistencies or model mismatches: {len(bad)}')
