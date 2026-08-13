"""Class groups of imaginary quadratic orders, as honest groups.

Until now every group-theoretic question about Cl(d) was answered by walking
horizontal isogeny cycles (qfs.qf_isog_cycle & co.).  This module supplies
the missing primitive -- Gauss/Dirichlet composition of binary quadratic
forms -- and wraps it in the AbGrp/AbGrElt paradigm of alg_classes.

Element tuples are the reduced primitive forms (a, b, c) of discriminant
d < 0, i.e. exactly the representatives qf_mod_gamma produces; constructing
an element from an arbitrary SL2(Z)-representative reduces it first.
Composition follows Cohen, A Course in Computational Algebraic Number
Theory, Algorithm 5.4.7 (plain composition; NUCOMP is unnecessary at our
discriminant sizes).  Non-maximal orders work throughout -- composition
never leaves the fixed discriminant d.

The prime bridge: class_from_prime(d, l) is the class of a form (l, b, *),
the same element x_l that generates the horizontal l-isogeny action.  A form
and its opposite represent the same primes, so a prime pins down a class
only up to inversion; primes_by_class therefore puts each prime in the
bucket of BOTH x_l and -x_l instead of choosing a sign.
"""

from functools import lru_cache

from nt import primesBetween, quad_rec, sqrt_mod_prime
from alg_classes import (AbGrp, AbGrElt, Subgroup, element_order,
                         invariant_factors, coset_decomposition,
                         quotient_group, abgrp_basis)
from qfs import (qf_mod_gamma, qf_disc, qf_is_prim, qf_in_fundom,
                 class_group_id, class_group_inv, get_qfs_strict)


                        ###############
                        # Composition #
                        ###############

def _xgcd(a: int, b: int) -> tuple[int, int, int]:
    """(g, u, v) with u*a + v*b = g = gcd(a, b) >= 0."""
    old_r, r = a, b
    old_u, u = 1, 0
    old_v, v = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_u, u = u, old_u - q * u
        old_v, v = v, old_v - q * v
    if old_r < 0:
        old_r, old_u, old_v = -old_r, -old_u, -old_v
    return old_r, old_u, old_v


def _compose_raw(qf1: tuple, qf2: tuple) -> tuple[int, int, int]:
    """Cohen Alg. 5.4.7 with no input validation (same-disc primitive forms
    assumed): the hot path for the group machinery, mirroring the
    _qf_reduce / qf_to_fun_dom split in qfs.py."""
    a1, b1, c1 = qf1
    a2, b2, c2 = qf2
    d = b1 * b1 - 4 * a1 * c1
    if a1 > a2:
        a1, b1, c1, a2, b2, c2 = a2, b2, c2, a1, b1, c1
    s = (b1 + b2) // 2
    n = b2 - s
    if a2 % a1 == 0:
        y1, dd = 0, a1
    else:
        dd, u, _ = _xgcd(a2, a1)
        y1 = u
    if s % dd == 0:
        y2, x2, d1 = -1, 0, dd
    else:
        d1, u, v = _xgcd(s, dd)
        x2, y2 = u, -v
    v1 = a1 // d1
    v2 = a2 // d1
    r = (y1 * y2 * n - x2 * c2) % v1
    a3 = v1 * v2
    b3 = b2 + 2 * v2 * r
    c3 = (b3 * b3 - d) // (4 * a3)
    return qf_mod_gamma((a3, b3, c3))


def qf_compose(qf1: tuple[int, int, int], qf2: tuple[int, int, int]) -> tuple[int, int, int]:
    """Gauss/Dirichlet composition of two primitive forms of one discriminant
    (Cohen, Algorithm 5.4.7), returned reduced via qf_mod_gamma."""
    a1, b1, c1 = qf1
    a2, b2, c2 = qf2
    if b1 * b1 - 4 * a1 * c1 != b2 * b2 - 4 * a2 * c2:
        raise ValueError('Forms must share a discriminant')
    if not (qf_is_prim(qf1) and qf_is_prim(qf2)):
        raise ValueError('Composition needs primitive forms')
    return _compose_raw(qf1, qf2)


                        ################
                        # Prime bridge #
                        ################

def class_from_prime(d: int, l: int):
    """The class x_l of a form (l, b, *) with b^2 = d (mod 4l): the class of
    a prime ideal above l, defined up to inversion.  Returns a reduced form,
    or None when l is inert or the form (l, b, *) is imprimitive (l divides
    the conductor)."""
    if d >= 0 or d % 4 > 1:
        raise ValueError(f'{d} is not a negative discriminant')
    if l == 2:
        r8 = d % 8
        if r8 == 5:
            return None                       # 2 inert
        b = {0: 0, 1: 1, 4: 2}[r8]
    else:
        if quad_rec(d, l) == -1:
            return None                       # l inert
        r = sqrt_mod_prime(d % l, l)
        b = r if (r - d) % 2 == 0 else l - r  # l odd: r, l-r have opposite parity
    c, rem = divmod(b * b - d, 4 * l)
    assert rem == 0, (d, l, b)
    qf = (l, b, c)
    if not qf_is_prim(qf):
        return None                           # l divides the conductor
    return qf_mod_gamma(qf)


def primes_by_class(d: int, N: int) -> dict:
    """{reduced class: sorted primes l <= N it represents}, every class of
    disc d present as a key.  Each prime lands in the bucket of x_l AND of
    -x_l: opposite forms represent the same primes, so this is the honest
    (sign-free) assignment."""
    out = {qf: [] for qf in class_group(d).elements_raw()}
    for l in primesBetween(2, N):
        x = class_from_prime(d, l)
        if x is None:
            continue
        out[x].append(l)
        xi = class_group_inv(x)
        if xi != x:
            out[xi].append(l)
    return out


def smallest_prime_by_class(d: int, cap: int = 1000, hard_cap: int = 1 << 20) -> dict:
    """{reduced class: smallest prime it represents}, growing the search cap
    until every class is hit (every primitive class represents infinitely
    many primes, so this terminates; hard_cap guards against typos)."""
    while True:
        buckets = primes_by_class(d, cap)
        if all(buckets.values()):
            return {qf: ls[0] for qf, ls in buckets.items()}
        if cap >= hard_cap:
            missing = [qf for qf, ls in buckets.items() if not ls]
            raise ValueError(f'No prime <= {cap} for classes {missing} of disc {d}')
        cap *= 4


                        ###############
                        # Group class #
                        ###############

class ClGrpElt(AbGrElt):
    """A class in ClassGroup(d).  Constructible from any integral form of
    discriminant d (any SL2(Z)-representative); stored reduced."""
    def __init__(self, qf, group: 'ClassGroup'):
        super().__init__(group.reduce(qf), group)

    def __eq__(self, other):
        if not isinstance(other, ClGrpElt) or other.grp is not self.grp:
            return NotImplemented
        return self.vec == other.vec

    def __hash__(self):
        return hash((id(self.grp), self.vec))

    @property
    def order(self) -> int:
        return element_order(self.grp, self.vec, self.grp.order)

    def primes(self, N: int) -> list[int]:
        """Primes <= N represented by this class (equivalently by -self)."""
        return [l for l in primesBetween(2, N)
                if class_from_prime(self.grp.disc, l) in (self.vec, (-self).vec)]


class ClassGroup(AbGrp):
    """Cl(d) for a negative discriminant d (non-maximal orders included).
    Element tuples are reduced primitive forms of discriminant exactly d.
    Prefer the cached factory class_group(d) so all elements of one
    discriminant share a single group object."""
    def __init__(self, d: int):
        if d >= 0 or d % 4 > 1:
            raise ValueError(f'{d} is not a negative discriminant')
        self.disc = d
        self._order = None
        self._elements_raw = None

        def membership(qf):
            return (isinstance(qf, tuple) and len(qf) == 3
                    and qf_in_fundom(qf) and qf_is_prim(qf)
                    and qf[1] * qf[1] - 4 * qf[0] * qf[2] == d)

        # addition_map is the unvalidated fast path (_compose_raw); the
        # public add_elements still membership-checks its inputs.
        super().__init__(membership, qf_mod_gamma(class_group_id(d)),
                         class_group_inv, _compose_raw)

    def __repr__(self):
        return f'ClassGroup({self.disc}), order {self.order}'

    def elements_raw(self) -> tuple:
        """The reduced primitive forms of disc d, enumerated once and cached."""
        if self._elements_raw is None:
            self._elements_raw = tuple(get_qfs_strict(self.disc))
        return self._elements_raw

    @property
    def order(self) -> int:
        """|Cl(d)|, counted from the reduced-forms enumeration (cheaper than
        the clgr_size_gen character sum; the two are cross-checked in
        scripts/validate_clgp.py)."""
        if self._order is None:
            self._order = len(self.elements_raw())
        return self._order

    def reduce(self, qf) -> tuple[int, int, int]:
        a, b, c = qf
        if b * b - 4 * a * c != self.disc:
            raise ValueError(f'{tuple(qf)} does not have discriminant {self.disc}')
        if not qf_is_prim((a, b, c)):
            raise ValueError(f'{tuple(qf)} is imprimitive: its class lives in '
                             'the class group of a smaller discriminant')
        return qf_mod_gamma((a, b, c))

    def element(self, qf) -> ClGrpElt:
        return ClGrpElt(tuple(qf), self)

    @property
    def identity(self) -> ClGrpElt:
        return ClGrpElt(self.zero_element, self)

    def elements(self) -> list[ClGrpElt]:
        return [ClGrpElt(qf, self) for qf in self.elements_raw()]

    def from_prime(self, l: int):
        """The element x_l (up to sign), or None if l is unusable."""
        qf = class_from_prime(self.disc, l)
        return None if qf is None else ClGrpElt(qf, self)

    # --- structure ---

    def basis(self) -> list[tuple[ClGrpElt, int]]:
        """Direct-product basis [(x_i, ord x_i)], orders non-increasing."""
        raw = abgrp_basis(self, self.elements_raw(), self.order)
        return [(ClGrpElt(t, self), o) for t, o in raw]

    def structure(self) -> tuple[int, ...]:
        """Invariant factors (d_1, ..., d_r), d_1 | ... | d_r."""
        return invariant_factors(o for _, o in self.basis())

    def subgroup(self, xs) -> Subgroup:
        """Subgroup generated by elements (ClGrpElt or raw form tuples)."""
        gens = [x.vec if isinstance(x, AbGrElt) else self.reduce(x) for x in xs]
        return Subgroup(self, gens)

    def cosets(self, H: Subgroup) -> list:
        return coset_decomposition(H, self.elements_raw())

    def quotient(self, H: Subgroup) -> AbGrp:
        return quotient_group(H, self.elements_raw())


@lru_cache(maxsize=None)
def class_group(d: int) -> ClassGroup:
    """Cached ClassGroup(d): one shared group object per discriminant."""
    return ClassGroup(d)


                #####################################
                # Prime-scored generating set search #
                #####################################

# A generating-set candidate is a +-pair {x, -x} tagged with the smallest
# prime representing it -- the natural currency of the isogeny side, where a
# generator costs the modular polynomial / Velu work of its prime.  The
# searches below look for direct-product bases among these candidates and
# rank them by a score on the primes used.

def default_score(primes_and_orders):
    """(largest prime, sum of primes): minimize the worst modular polynomial
    first, total prime weight second."""
    ls = [l for l, _ in primes_and_orders]
    return (max(ls), sum(ls)) if ls else (0, 0)


def _prime_pairs_lazy(G: ClassGroup, bound: int, known: dict, scanned: int) -> int:
    """Extend known ({pair key: smallest prime}) with the usable primes in
    (scanned, bound]; returns the new scanned bound.  Pool built in increasing
    prime order, so known always holds ALL pairs whose smallest prime is <=
    the scanned bound."""
    for l in primesBetween(scanned + 1, bound):
        x = class_from_prime(G.disc, l)
        if x is None or x == G.zero_element:
            continue
        key = min(x, class_group_inv(x))
        if key not in known:
            known[key] = l
    return bound


def _first_prime_in_classes(d: int, targets, skip, cap: int = 1000,
                            hard_cap: int = 1 << 20):
    """Smallest prime (not in skip) whose class or inverse class lies in
    targets, or None past hard_cap."""
    lo = 2
    while lo <= hard_cap:
        for l in primesBetween(lo, cap):
            if l in skip:
                continue
            x = class_from_prime(d, l)
            if x is not None and (x in targets or class_group_inv(x) in targets):
                return l
        lo, cap = cap + 1, cap * 4
    return None


def prime_generator_candidates(G: ClassGroup, prime_cap: int = None) -> list:
    """[(l, qf)] with qf the reduced class of the smallest prime l
    representing the +-pair {x, -x}, sorted by l.  Grows the prime bound
    lazily until the collected pairs span Cl(d) (so it stays cheap when small
    primes already generate); prime_cap forces a fixed bound instead."""
    known: dict = {}
    if prime_cap is not None:
        _prime_pairs_lazy(G, prime_cap, known, 0)
    else:
        scanned, bound = 0, 200
        while True:
            scanned = _prime_pairs_lazy(G, bound, known, scanned)
            if len(Subgroup(G, list(known))) == G.order:
                break
            if bound > (1 << 20):
                raise ValueError(f'primes to {bound} do not span Cl({G.disc})')
            bound *= 4
    return sorted(((l, key) for key, l in known.items()), key=lambda t: t[0])


def _prime_pairs_from_pool(G: ClassGroup, pool) -> dict:
    """{pair key: smallest prime in pool representing it}."""
    known: dict = {}
    for l in sorted(pool):
        x = class_from_prime(G.disc, l)
        if x is None or x == G.zero_element:
            continue
        key = min(x, class_group_inv(x))
        if key not in known:
            known[key] = l
    return known


def prime_scored_bases(G: ClassGroup, prime_cap: int = None, score=None,
                       top: int = 1, pool=None) -> list:
    """The top-scoring direct-product bases drawn from prime-represented
    elements.  Returns [(score, [(l, qf, ord), ...])], best first.

    score maps a list of (prime, order) pairs to a comparable; it must be
    monotone under extending the list (the default is), since partial bases
    are pruned against the current top scores.  Candidate primes: an explicit
    pool if given (e.g. modular_prime_pool()), else all primes to prime_cap,
    else a lazily growing bound -- in which case the result is optimal for
    any score dominated by the largest prime used."""
    if score is None:
        score = default_score
    n = G.order

    def search(known: dict) -> list:
        found = []
        # a spanning pool need not contain a basis ({(1,1),(1,2)} spans
        # Z2 x Z4 basis-free), so DFS decides; spanning is just a cheap gate
        if n > 1 and len(Subgroup(G, list(known))) != n:
            return found
        cands = sorted(((l, key) for key, l in known.items()),
                       key=lambda t: t[0])
        orders = {qf: element_order(G, qf, n) for _, qf in cands}

        def worse_than_top(sc):
            return len(found) == top and sc >= found[-1][0]

        def extend(basis, span, start):
            size = len(span)
            if size == n:
                sc = score([(l, o) for l, _, o in basis])
                if not worse_than_top(sc):
                    found.append((sc, list(basis)))
                    found.sort(key=lambda t: t[0])
                    del found[top:]
                return
            for i in range(start, len(cands)):
                l, qf = cands[i]
                o = orders[qf]
                if (size * o > n) or (n % (size * o)) or (qf in span):
                    continue
                if worse_than_top(score([(l2, o2) for l2, _, o2 in basis]
                                        + [(l, o)])):
                    continue
                mults, kt, direct = [], qf, True
                for _ in range(o - 1):
                    if kt in span:
                        direct = False
                        break
                    mults.append(kt)
                    kt = G.addition_map(kt, qf)
                if not direct:
                    continue
                new_span = set(span)
                for kt in mults:
                    new_span.update(G.addition_map(s, kt) for s in span)
                extend(basis + [(l, qf, o)], new_span, i + 1)

        extend([], {G.zero_element}, 0)
        return found

    if pool is not None:
        return search(_prime_pairs_from_pool(G, pool))
    known: dict = {}
    scanned, bound = 0, (200 if prime_cap is None else prime_cap)
    while True:
        scanned = _prime_pairs_lazy(G, bound, known, scanned)
        found = search(known)
        if found or prime_cap is not None or bound > (1 << 20):
            return found
        bound *= 4


def rigid_lset(G: ClassGroup, prime_cap: int = None, score=None, pool=None):
    """A rigid l-set for Cl(d) computed with the group law: a best-scoring
    prime basis, plus (when there are >= 2 generators of order > 2) a pinning
    prime representing a signed sum of the long generators.

    Returns {'ls_rig': (longs..., l_sum, twos...), 'basis': ..., 'pin': ...}
    matching the tuple contract of ecqf_bij / compute_bijection_zn, or None
    if no basis or pin is available (only possible with an explicit
    prime_cap or pool).  Pass pool=modular_prime_pool() to search under the
    j-side constraint like disc_rigid_lset_search does.  This is the fast
    path of that search's simple case; sib/lift/power descriptors and
    relaxed pins stay over there for now."""
    bases = prime_scored_bases(G, prime_cap, score=score, top=1, pool=pool)
    if not bases:
        return None
    _, basis = bases[0]
    longs = [(l, qf, o) for l, qf, o in basis if o > 2]
    twos = [(l, qf, o) for l, qf, o in basis if o == 2]
    pin = None
    if len(longs) >= 2:
        sums = {G.zero_element}
        for _, qf, _o in longs:
            sums = ({G.addition_map(s, qf) for s in sums}
                    | {G.addition_map(s, class_group_inv(qf)) for s in sums})
        sums.discard(G.zero_element)
        used = {l for l, _, _ in basis}
        if pool is not None:
            pin = next((l for l in sorted(pool) if l not in used
                        and (x := class_from_prime(G.disc, l)) is not None
                        and (x in sums or class_group_inv(x) in sums)), None)
        else:
            pin = _first_prime_in_classes(G.disc, sums, used,
                                          cap=prime_cap or 1000,
                                          hard_cap=prime_cap or (1 << 20))
        if pin is None:
            return None
    ls_rig = (tuple(l for l, _, _ in longs)
              + ((pin,) if pin is not None else ())
              + tuple(l for l, _, _ in twos))
    return {'ls_rig': ls_rig, 'basis': basis, 'pin': pin}
