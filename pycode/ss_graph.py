# The true supersingular l-isogeny graph over F_{p^2}.
#
# Vertices are ALL supersingular j-invariants in F_{p^2} (the F_p-rational
# ones used by the volcano pages plus the conjugate pairs); the directed edge
# multiplicity adj[i][k] is the multiplicity of j_k as a root of
# Phi_l(j_i, Y), which equals the number of order-l subgroups C of E_{j_i}
# with j(E/C) = j_k, so every row sums to l + 1.
#
# Structure carried along with the graph:
#   - rationality of each vertex (j in F_p);
#   - the Frobenius involution j -> j^p on vertices and on undirected edges
#     (the number of conjugate vertex pairs equals the genus of X_0(p)^+);
#   - the Eichler weighted symmetry  e_k * adj[i][k] == e_i * adj[k][i]
#     (dualizing is a bijection on isogenies and #Isog(E_i, E_k) =
#     adj[i][k] * |Aut(E_k)|), e_j = |Aut(E_j)| / 2: 2 at j = 1728, 3 at
#     j = 0, else 1 — the only failure of plain symmetry.
#
# F_{p^2} elements are int pairs (u, v) meaning u + v*sqrt(s), s the smallest
# nonresidue mod p — the same convention as pell_duality.py, whose F_q
# toolkit this module reuses.  Kept synchronous and modest: fine for p up to
# a few hundred and the l with stored classical modular polynomials.

import json
from functools import lru_cache
from pathlib import Path

import numpy as np

from pell_duality import (
    smallest_nonresidue, ss_vertex_set, phi_fiber_poly,
    qeval, qdivmod, qnorm, qdeg, fq_conj,
)

_DATA_DIR = Path(__file__).parent / 'data'


@lru_cache(maxsize=1)
def available_ls():
    """Primes l with a stored classical modular polynomial Phi_l."""
    with open(_DATA_DIR / 'classical_modpolys.json') as f:
        raw = json.load(f)
    return tuple(sorted(int(k) for k in raw))


def jstr(j, s):
    """Display string for j = (u, v) in F_{p^2} = F_p(sqrt(s))."""
    u, v = j
    if v == 0:
        return str(u)
    return f'{u}+{v}√{s}'


def aut_weight(j, p):
    """e_j = |Aut(E_j)| / 2 for supersingular j (p >= 5)."""
    if j == (0, 0):
        return 3
    if j == (1728 % p, 0):
        return 2
    return 1


def _linear(root, p):
    """The monic linear q-poly Y - root."""
    return qnorm(np.array([[(-root[0]) % p, 1],
                           [(-root[1]) % p, 0]], dtype=np.int64), p)


def ss_isogeny_graph(p, l):
    """The supersingular l-isogeny graph over F_{p^2}, with multiplicities.

    Returns a dict:
      p, l, s     : the prime, the isogeny degree, the nonresidue defining F_q
      vertices    : sorted list of all supersingular j in F_{p^2}, as (u, v)
      rational    : rational[i] is True iff vertices[i] lies in F_p
      frob        : frob[i] = index of the Frobenius conjugate of vertices[i]
      adj         : adj[i][k] = # order-l subgroups of E_{j_i} landing on j_k
                    (= multiplicity of j_k in Phi_l(j_i, Y)); rows sum to l+1
      edges       : undirected edges, one dict per pair i <= k with a nonzero
                    multiplicity: {i, k, m (i->k), m_rev (k->i), frob_edge
                    (index of the conjugate edge), frob_class ('fixed' if the
                    edge equals its conjugate, else 'swapped')}
    """
    if p < 5:
        raise ValueError('p must be a prime >= 5')
    if p == l:
        raise ValueError('p and l must be distinct')
    s = smallest_nonresidue(p)
    vertices = ss_vertex_set(p, s)
    n = len(vertices)
    idx = {j: i for i, j in enumerate(vertices)}

    # Directed multiplicities: peel the roots of Phi_l(j_i, Y) off the ss
    # vertex set (supersingular is isogeny-invariant, so the fiber splits
    # completely there).
    adj = [[0] * n for _ in range(n)]
    for i, j in enumerate(vertices):
        remaining = phi_fiber_poly(l, j, p, s)
        assert qdeg(remaining) == l + 1, (p, l, j, qdeg(remaining))
        for k, jk in enumerate(vertices):
            while qdeg(remaining) > 0 and qeval(remaining, jk, p, s) == (0, 0):
                quo, rem = qdivmod(remaining, _linear(jk, p), p, s)
                assert rem.shape[1] == 0
                remaining = quo
                adj[i][k] += 1
        assert qdeg(remaining) == 0, (p, l, j, 'fiber did not split on ss set')
        assert sum(adj[i]) == l + 1

    rational = [j[1] == 0 for j in vertices]
    frob = [idx[fq_conj(j, p)] for j in vertices]

    # Structural checks: Frobenius equivariance and Eichler weighted symmetry.
    e = [aut_weight(j, p) for j in vertices]
    for i in range(n):
        for k in range(n):
            assert adj[i][k] == adj[frob[i]][frob[k]], 'not Frobenius-equivariant'
            assert e[k] * adj[i][k] == e[i] * adj[k][i], 'Eichler symmetry fails'

    edges = []
    edge_idx = {}
    for i in range(n):
        for k in range(i, n):
            if adj[i][k]:
                edge_idx[(i, k)] = len(edges)
                edges.append({'i': i, 'k': k,
                              'm': adj[i][k], 'm_rev': adj[k][i]})
    for eid, ed in enumerate(edges):
        key = tuple(sorted((frob[ed['i']], frob[ed['k']])))
        ed['frob_edge'] = edge_idx[key]
        ed['frob_class'] = 'fixed' if ed['frob_edge'] == eid else 'swapped'

    return {'p': p, 'l': l, 's': s, 'vertices': vertices,
            'rational': rational, 'frob': frob, 'adj': adj, 'edges': edges}


def kernel_selfcheck(p, l, seed=0):
    """Recompute the adjacency at kernel level and compare.

    Enumerates all l+1 kernel polynomials at every vertex (pell_duality's
    Velu machinery) and counts codomain j-invariants; must equal the
    Phi-multiplicity adjacency exactly."""
    from pell_duality import (std_model, kernels_at_vertex, velu_codomain,
                              fq_j_invariant)
    g = ss_isogeny_graph(p, l)
    s, vertices, n = g['s'], g['vertices'], len(g['vertices'])
    idx = {j: i for i, j in enumerate(vertices)}
    for i, j in enumerate(vertices):
        a, b = std_model(j, p, s)
        counts = [0] * n
        for K in kernels_at_vertex(a, b, l, p, s, seed=seed):
            av, bv = velu_codomain(a, b, K, l, p, s)
            counts[idx[fq_j_invariant(av, bv, p, s)]] += 1
        assert counts == g['adj'][i], (p, l, j, counts, g['adj'][i])
    return True


if __name__ == '__main__':
    import time
    for p, l in [(11, 2), (13, 2), (23, 3), (31, 5), (37, 2), (47, 3)]:
        t0 = time.time()
        kernel_selfcheck(p, l)
        g = ss_isogeny_graph(p, l)
        n = len(g['vertices'])
        pairs = sum(1 for i in range(n) if i < g['frob'][i])
        print(f'p={p:3d} l={l:2d}: {n:2d} vertices '
              f'({sum(g["rational"])} rational, {pairs} pairs), '
              f'{len(g["edges"])} edges — kernel selfcheck OK '
              f'({time.time() - t0:.1f}s)')
