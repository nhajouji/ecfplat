from misctools import *

####################
# Adjacency Matrix #
####################

def nbrdata_to_isomat(nbrdata:dict,verts_ordered:list):
    v_to_i = {v:i for i,v in enumerate(verts_ordered)}
    assert len({v for v in v_to_i if v not in nbrdata}) == 0
    dim = len(v_to_i)
    mat = [[0 for _ in range(dim)] for _ in range(dim)]
    for i1,v1 in enumerate(verts_ordered):
        for v2 in nbrdata[v1]:
            i2 = v_to_i[v2]
            mat[i1][i2]+=1
    return mat

##########################
# Tree search algorithms #
##########################


def nbrdata_tree_search(nbrdata:dict,l_gens:tuple[int],leaf_cands:list,tree_root):
    chains = [[tree_root]]
    for l in l_gens:
        if l not in nbrdata:
            raise ValueError(f'Missing data for l = {l}')
        lnbrs = nbrdata[l]
        newchains = []
        for chain in chains:
            leaf = chain[-1]
            if leaf not in lnbrs or len(lnbrs[leaf])==0:
                raise ValueError(f'Isogeny data missing for {leaf} for deg {l}')
            for dsc in lnbrs[leaf]:
                if dsc not in chain:
                    newchains.append(chain+[dsc])
        chains = newchains
    chains_admiss = [ch for ch in chains if ch[-1] in leaf_cands]
    return [{l_gens[:i]:v for i,v in enumerate(chain)} for chain in chains_admiss]

def nbrdata_tree_search_XCyc(nbrdata:dict,l_gens:tuple[int],leaf_cands:list,tree_root):
    # TODO (known limitation, fine for now): the radical-word search re-walks
    # the tree over CYCLIC ROTATIONS of l_gens only.  Every subset of <= 3
    # generators is a prefix of some rotation, so the contract (determine F at
    # ALL radical words) holds for up to 3 generators of order > 2; with 4+,
    # non-adjacent supports like {b1, b3} are prefixes of no rotation and
    # those radical words are never determined.  Fix: walk one ordering per
    # radical-word support (any ordering whose prefix is the support).  Not
    # urgent: no discriminant with |d| < 64000 has 4 such generators (first
    # 3-generator cases appear near p ~ 16000).
    chains = [[tree_root]]
    l_gens_original_order = {l:i for i,l in enumerate(l_gens)}
    for l in l_gens:
        if l not in nbrdata:
            raise ValueError(f'Missing data for l = {l}')
        lnbrs = nbrdata[l]
        newchains = []
        for chain in chains:
            leaf = chain[-1]
            if leaf not in lnbrs or len(lnbrs[leaf])==0:
                raise ValueError(f'Isogeny data missing for {leaf} for deg {l}')
            for dsc in lnbrs[leaf]:
                if dsc not in chain:
                    newchains.append(chain+[dsc])
        chains = newchains
    chains_admiss = [ch for ch in chains if ch[-1] in leaf_cands]
    if len(chains_admiss)== 0:
        raise ValueError('No chains found')
    ch0 = chains_admiss[0]
    leaf0 = ch0[-1]
    return [{sort_tuple(v,l_gens_original_order):x for v,x in nbrdata_tree_search(nbrdata,rotate_tuple(l_gens,k),[leaf0],tree_root)[0].items()}
            for k in range(len(l_gens))]


def nbrdata_tree_search_XC1(nbrdata:dict,l_gens:tuple[int],leaf_cands:list,tree_root):
    return merge_dicts(nbrdata_tree_search_XCyc(nbrdata=nbrdata,l_gens=l_gens,leaf_cands=leaf_cands,tree_root=tree_root))

def nbrdata_tree_search_zn(nbrdata:dict,l_gens:tuple[int],leaf_cands:list,tree_root):
    orig_d = nbrdata_tree_search_XC1(nbrdata = nbrdata,
                                     l_gens=l_gens,
                                     leaf_cands=leaf_cands,
                                     tree_root=tree_root)
    return {tuple([int(l in lt) for l in l_gens]):orig_d[lt] for lt in orig_d}


"""
Labeling Cayley graphs of finite abelian groups.

Setup
-----
A = Z_{n_1} x ... x Z_{n_r} with independent generators e_1, ..., e_r,
each of order n_i > 2.  Group elements are represented by integer tuples
(x_1, ..., x_r) with 0 <= x_i < n_i.

We are given:
  * nbrdata: a list [D_1, ..., D_r] of dictionaries.  D_i maps each v in V
    to the unordered pair [F(x - e_i), F(x + e_i)] where v = F(x).
  * initF: a dictionary giving F on the unit cube, i.e. with keys all
    0/1-tuples of length r (so 2^r keys) and values in V.

extend_bijection(initF, nbrdata) returns the unique extension of initF to a
bijection F : A -> V compatible with all the D_i, as a dictionary whose keys
are tuples (x_1, ..., x_r).

Convention: coordinate i of a key (0-indexed) corresponds to direction
nbrdata[i].  Note this differs from oriented_l12page_nbd in the original
notebook, where the output key (i, j) had i indexing the *second* direction.

The orders n_i are not needed as input: they are discovered as the cycle
lengths, and reported by group_orders().
"""


def oriented_cycle_from_nbdata(v01, nbrdata):
    """Walk the direction-i cycle starting v0 -> v1.

    Returns the list [F(x), F(x+e_i), F(x+2e_i), ..., F(x+(n_i-1)e_i)]
    where v0 = F(x), v1 = F(x+e_i).  Linear time: at each step the next
    vertex is the neighbor of the current one that is not the previous one
    (well-defined because each vertex has exactly two distinct neighbors,
    guaranteed by ord(e_i) > 2).
    """
    v0, v1 = v01
    if v0 not in nbrdata or v1 not in nbrdata:
        raise KeyError("Vertices must be in dictionary")
    if v1 not in nbrdata[v0]:
        raise ValueError("Vertices must be neighbors")
    cyc = [v0, v1]
    prev, cur = v0, v1
    while True:
        a, b = nbrdata[cur]
        nxt = b if a == prev else a
        if nxt == v0:
            return cyc
        cyc.append(nxt)
        prev, cur = cur, nxt

def cycle_from_neighbor_data(v0,nbrdata):
    if v0 not in nbrdata:
        raise KeyError('Vertices must be in dictionary')
    v0_nbrs = nbrdata[v0]
    if len(v0_nbrs)<=1:
        return [v0]+v0_nbrs
    else:
        v1 = v0_nbrs[0]
        return oriented_cycle_from_nbdata(v01=(v0,v1),nbrdata=nbrdata)

def extend_bijection_zn1(initF, nbrdata):
    """Extend initF (defined on the 0/1 cube) to all of A.

    Induction on r:
      r = 1: a single oriented cycle F(0) -> F(e_1).
      r > 1: split the cube into its two faces x_r = 0 and x_r = 1, fill
             each (r-1)-dimensional slab recursively, then for every point
             x of the slab x_r = 0 run the cycle F(x) -> F(x + e_r).
    """
    r = len(nbrdata)
    _check_cube(initF, r)
    return _fill(initF, nbrdata)


def _fill(cube, nbrs):
    r = len(nbrs)
    if r == 1:
        cyc = oriented_cycle_from_nbdata((cube[(0,)], cube[(1,)]), nbrs[0])
        return {(m,): v for m, v in enumerate(cyc)}
    face0 = {k[:-1]: v for k, v in cube.items() if k[-1] == 0}
    face1 = {k[:-1]: v for k, v in cube.items() if k[-1] == 1}
    slab0 = _fill(face0, nbrs[:-1])
    slab1 = _fill(face1, nbrs[:-1])
    F = {}
    n_last = None
    for key, v0 in slab0.items():
        cyc = oriented_cycle_from_nbdata((v0, slab1[key]), nbrs[-1])
        if n_last is None:
            n_last = len(cyc)
        elif len(cyc) != n_last:
            raise ValueError(
                f"Inconsistent cycle lengths in direction {r - 1}: "
                f"{n_last} vs {len(cyc)}"
            )
        for m, v in enumerate(cyc):
            F[key + (m,)] = v
    return F


def _check_cube(initF, r):
    from itertools import product
    expected = set(product((0, 1), repeat=r))
    got = set(initF)
    if got != expected:
        missing = expected - got
        extra = got - expected
        raise ValueError(
            f"initF must have exactly the 0/1 tuples of length {r} as keys; "
            f"missing {sorted(missing)}, unexpected {sorted(extra)}"
        )


def group_orders(F):
    """Recover (n_1, ..., n_r) from a completed labeling."""
    r = len(next(iter(F)))
    return tuple(1 + max(k[i] for k in F) for i in range(r))


def verify_bijection(F, nbrdata):
    """Sanity check: F is injective, and for every x and every direction i,
    {F(x + e_i), F(x - e_i)} equals the pair nbrdata[i][F(x)].
    Raises AssertionError on failure."""
    vals = list(F.values())
    assert len(set(vals)) == len(vals), "F is not injective"
    ns = group_orders(F)
    for x, v in F.items():
        for i, D in enumerate(nbrdata):
            xp = x[:i] + ((x[i] + 1) % ns[i],) + x[i + 1:]
            xm = x[:i] + ((x[i] - 1) % ns[i],) + x[i + 1:]
            assert {F[xp], F[xm]} == set(D[v]), (
                f"Neighbor mismatch at x={x}, direction {i}"
            )
    return True

##########################################
# Relaxed pin: powered cube -> unit cube #
##########################################
# The pinning element of a rigid l-set may be relaxed from the plain sum of
# the long generators to  sum(c_i * b_i) + tau,  with full support,
# 0 < c_i < ord(b_i)/2, and tau a word in the order-2 generators.  The
# assembler then walks the isogeny tree with the c_i-step graphs (leaf set
# {sum +-c_i b_i}, collision-free by c_i < ord/2), after stripping tau from
# the pin's neighbours of the root (order-2 steps are orientation-free, and
# -tau = tau, so both signs strip correctly).  The resulting 0/1 cube has
# 1-direction i meaning +c_i b_i; cube_from_powered converts it to the unit
# cube so extend_bijection_zn1 applies unchanged.

def _is_pin_label(x):
    """Structural check for a relaxed-pin label ('pin', desc, pairs, taus)."""
    return isinstance(x, tuple) and len(x) == 4 and x[0] == 'pin'

def _offset_first_step(nbrs, v0, vc, c):
    """v0 and vc sit at distance c (0 < c < n/2) on the direction cycle nbrs:
    return the 1-step neighbour of v0 on the vc side.  The bound c < n/2 is
    what makes the side well-defined."""
    cyc = cycle_from_neighbor_data(v0, nbrs)
    n = len(cyc)
    if cyc[c % n] == vc:
        return cyc[1 % n]
    assert cyc[(-c) % n] == vc, 'offset vertex not at the stated distance'
    return cyc[-1]

def cube_from_powered(cube_pow, base_nbrs, cs):
    """Convert a 0/1 cube whose 1 in direction i means +c_i e_i into the unit
    cube (+e_i), fixing one coordinate at a time.  base_nbrs[i] is the 1-step
    graph of direction i.  After fixing coordinate i, the pair (k0, k1) along
    any later coordinate j still differs by c_j e_j, so the induction holds."""
    cube = dict(cube_pow)
    for i, c in enumerate(cs):
        if c == 1:
            continue
        new = {}
        for k, v in cube.items():
            if k[i] == 0:
                new[k] = v
            else:
                v0 = cube[k[:i] + (0,) + k[i + 1:]]
                new[k] = _offset_first_step(base_nbrs[i], v0, v, c)
        cube = new
    return cube


#######################
# 2-torsion extension #
#######################

def extend_bijection_zn2_x1(bij_initial:dict,l2neighbordata:dict):
    newdict = {}
    for t in bij_initial:
        vt = bij_initial[t]
        vt_nbrs = [v1 for v1 in l2neighbordata[vt] if v1 not in bij_initial]
        assert len(vt_nbrs)==1
        vt1 = vt_nbrs[0]
        newdict[ext_tup(t,0)]=vt
        newdict[ext_tup(t,1)] = vt1
    return newdict

def extend_bijection_zn2(bij0:dict,neighbor_data_all:dict,ls_t2:tuple[int]):
    for l in ls_t2:
        bij0 = extend_bijection_zn2_x1(bij0,neighbor_data_all[l])
    return bij0

#########################
# End to end bijection  #
#########################

def l2_split(neighbordata,ls,root):
    l2s = []
    lbigs = []
    for l in ls:
        if len(neighbordata[l][root])==1:
            l2s.append(l)
        elif len(neighbordata[l][root])==2:
            lbigs.append(l)
    return tuple(lbigs),tuple(l2s)
    
    
def compute_bijection_zn(ls_rig:tuple[int],neighbordata:dict,root):
    if len(ls_rig)==0:
        return {():root}
    # We split the set of l's into two groups - ls1 contains graph labels where cycles
    # have length greater than 2, and ls2 contains graph labels where cycles have length 2
    ls1,ls2 = l2_split(neighbordata=neighbordata,ls=ls_rig,root=root)
    # We just need to take care of the long cycles and use previous func to ext
    if len(ls1)==0:
        # There are no long cycles so we just extend the trivial dictionary 
        bij={():root}
        return extend_bijection_zn2(bij,neighbordata,ls2)
    elif len(ls1) == 1:
        # If theres only one long cycle, we don't have to worry about orientation
        # We compute the long cycle first and pick the orientation randomly
        # Then we extend using short cycles to complete th ebijection
        l = ls1[0]
        nbrdata_l = neighbordata[l]
        v0 = root
        assert root in nbrdata_l and len(nbrdata_l[root]) ==2
        v1 = nbrdata_l[root][0]
        bij = extend_bijection_zn1({(0,):v0,(1,):v1},[nbrdata_l])
        return extend_bijection_zn2(bij,neighbordata,ls2)
    else:
        pin = ls1[-1]
        leaf_cands = neighbordata[pin][root]
        if _is_pin_label(pin):
            # relaxed pin ('pin', pdesc, pairs, taus): pairs[i] = (key of the
            # c_i-step graph of long generator i, c_i), taus = order-2 labels
            # to strip from the pin's neighbours before the tree match.
            _, pdesc, pairs, taus = pin
            assert len(pairs) == len(ls1) - 1, 'pin pairs misaligned with longs'
            for t in taus:
                leaf_cands = [neighbordata[t][v][0] for v in leaf_cands]
            leaf_cands = list(dict.fromkeys(leaf_cands))
            cubedata = nbrdata_tree_search_zn(nbrdata=neighbordata,
                                              l_gens=tuple(pk for pk, c in pairs),
                                              leaf_cands=leaf_cands,
                                              tree_root=root)
            cubedata = cube_from_powered(cubedata,
                                         [neighbordata[l] for l in ls1[:-1]],
                                         [c for pk, c in pairs])
            bij = extend_bijection_zn1(cubedata, [neighbordata[l] for l in ls1[:-1]])
            return extend_bijection_zn2(bij, neighbordata, ls2)
        cubedata = nbrdata_tree_search_zn(nbrdata=neighbordata,
                                           l_gens=ls1[:-1],
                                           leaf_cands=leaf_cands,
                                           tree_root=root)
        bij = extend_bijection_zn1(cubedata,[neighbordata[l] for l in ls1[:-1]])
        return extend_bijection_zn2(bij,neighbordata,ls2)
