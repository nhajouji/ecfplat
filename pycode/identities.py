"""Class-number identities and mass counts.

The twisted Euler phi (the conductor factor in h(f^2 d) / h(d)), class-number
sums over conductor divisors (clgr_size_gen / clgr_sum, per the Hurwitz class
number), the supersingular trace count supsingtrace (Eichler mass formula
territory -- the 2p+2 identity of the p-l duality work), discriminant
closures, and point counts x0l_fp_card for X_0(l)(F_p).

Convention: everything here counts by pure character/divisor sums with no
form enumeration, so it doubles as an independent cross-check of qfs/clgp
(and is what the future test suite checks them against).
"""

from nt import gcd,quad_rec,gen_quad_symb,primefact,divisors,discfac,no_odd_prime_facs

def twisted_phi(d:int,m:int)->int:
    if m < 0:
        return 0
    pfm = primefact(m)
    phim = 1
    for p in pfm:
        phim*=(p-quad_rec(d,p))*(p**(pfm[p]-1))
    if d == -3 and m > 1:
        return phim//3
    elif d == - 4 and m > 1:
        return phim//2
    else:
        return phim

def twisted_phi_sum(d:int,m:int)->int:
    divsm = divisors(m)
    return sum([twisted_phi(d,m) for m in divsm])

# This computes the size of the class group of a fundamental discriminant d
# Mainly used as a helper function for clgr_size_gen
def clgr_size_fd(d:int)->int:
    if abs(d)<13:
        return 1
    den = 2-quad_rec(d,2)
    return sum([gen_quad_symb(d,l) for l in range(1,(abs(d)+1)//2) if gcd(l,d)==1])//den

# This returns the size of the class group of a (possibly nonmaximal) complex order
# from its discriminant. Equivalent to computing the size of get_qfs_strict
def clgr_size_gen(dc:int)->int:
    if abs(dc)<13:
        return 1
    d,c = discfac(dc)
    return clgr_size_fd(d)*twisted_phi(d,c)

# This returns the sum of the class group sizes of all complex orders that contains
# the order of discriminant d.
# Equivalent to computing the size of get_qfs_all
def clgr_sum(dc:int):
    d,c = discfac(dc)
    return sum([clgr_size_gen(d*c0*c0) for c0 in divisors(c)])

# This gives the (adjusted)
# exponent of the maximal 2-group in the class group of the complex order
# of discriminant d
def clgr_2len(d:int):
    if d % 4 == 1:
        return no_odd_prime_facs(d)
    n = d//(-4)
    r = no_odd_prime_facs(d)
    if n % 4 == 3:
        return r
    elif n % 4 in [1,2] or n % 8 == 4:
        return r+1
    else:
        return r+2

    
def clgr2_size(d):
    return 2**(clgr_2len(d)-1)

###############################
# Modular curves and SS trace #

def supsingtrace(p,l):
    trace = 0
    d3seen = False
    d4seen = False
    a = 0
    while a*a < 4*l:
        d,c = discfac(a*a-4*l)
        qr = quad_rec(d,p)
        if qr < 1:
            while c % p == 0:
                c = c // p
            d0 = d*c*c
            h = clgr_sum(d0)
            if d % p == 0 or d % l == 0:
                if d == -3: 
                    if not d3seen:
                        trace += h
                        d3seen = True
                    else:
                        trace+=h-1
                elif d == -4:
                    if not d4seen:
                        trace += h
                        d4seen = True
                    else:
                        trace+=h-1
                else:  
                    trace += h
            else:
                if d == -3: 
                    if not d3seen:
                        trace += 2*h
                        d3seen = True
                    else:
                        trace+=2*(h-1)
                elif d == -4:
                    if not d4seen:
                        trace += 2*h
                        d4seen = True
                    else:
                        trace+=2*(h-1)
                else:  
                    trace += 2*h
        a+=1                
    return trace

def disc_closure(d:int):
    if d >= 0 or d % 4 >1:
        raise ValueError(f'{d} must be a negative discriminant')
    d0, c = discfac(d)
    return [d0*(c0**2) for c0 in divisors(c)]

def ecfp_all_endo_discs(p:int):
    ds = []
    a = 0
    while a * a < 4*p:
        ds+=disc_closure(a*a-4*p)
        a+=1
    return ds

def ecfp_disc_l_isocts(p:int,l:int):
    ds_all = ecfp_all_endo_discs(p)
    cntrs = {d:(1+quad_rec(d,l)) for d in ds_all}
    for d in cntrs:
        if d*(l**2) in cntrs:
            if d == -3:
                cntrs[d]+=(l-quad_rec(d,l))//3
            elif d == - 4:
                cntrs[d]+=(l-quad_rec(d,l))//2
            else:
                cntrs[d]=l+1
    return cntrs

def x0l_fp_card(p:int,l:int):
    # Initialize ct at 2 for the pair of nodes
    ct = 2
    # Collect data for each isogeny class
    data = ecfp_disc_l_isocts(p,l)
    for d in data:
        if d %p != 0:
            ct+=data[d]*clgr_size_gen(d)
    if quad_rec(-p,l)==1:
        ct+=clgr_sum(-4*p)
    return ct

#############################################
# Genera of X_0(N), X_0(N)^+ and X_0(N)^*   #
#############################################
# References: Shimura Prop. 1.40/1.43 (genus of X_0(N)); Ogg, "Hyperelliptic
# modular curves" (Bull. SMF 1974) and Kluit for the fixed points of the
# Atkin-Lehner involutions w_Q; Riemann-Hurwitz for the quotients.

def _mu_index(N:int)->int:
    """[SL_2(Z) : Gamma_0(N)] = N prod_{p|N} (1 + 1/p)."""
    mu = N
    for q in primefact(N):
        mu = mu // q * (q + 1)
    return mu

def _nu2(N:int)->int:
    """Number of elliptic points of order 2 on X_0(N)."""
    if N % 4 == 0:
        return 0
    n = 1
    for q in primefact(N):
        n *= 1 + quad_rec(-4, q)     # (-1/q) via Kronecker (-4/q)
    return n

def _nu3(N:int)->int:
    """Number of elliptic points of order 3 on X_0(N)."""
    if N % 9 == 0:
        return 0
    n = 1
    for q in primefact(N):
        n *= 1 + quad_rec(-3, q)
    return n

def _nuinf(N:int)->int:
    """Number of cusps of X_0(N) = sum_{d | N} phi(gcd(d, N/d))."""
    tot = 0
    for d in divisors(N):
        g = gcd(d, N // d)
        tot += twisted_phi(1, g) if g > 1 else 1   # Euler phi via twisted_phi(1, .)
    return tot

def genus_x0(N:int)->int:
    """Genus of the modular curve X_0(N)."""
    if N == 1:
        return 0
    g12 = 12 + _mu_index(N) - 3 * _nu2(N) - 4 * _nu3(N) - 6 * _nuinf(N)
    assert g12 % 12 == 0, N
    return g12 // 12

def _hw(d:int)->int:
    """Class number of the order of discriminant d, with h(-3)=h(-4)=1
    (unweighted; the weights are put in by hand where they matter)."""
    return clgr_size_gen(d)

def al_fixed_points(Q:int, N:int)->int:
    """Number of fixed points of the Atkin-Lehner involution w_Q on X_0(N),
    for Q an exact divisor of N (Q || N), following Fricke / Ogg (Bull. SMF
    1974, Prop. 3) / Kluit.  Fixed points are CM points of discriminant -4Q,
    and also -Q when Q = 3 mod 4 (and -4, -8 for Q = 2; -3, -12 for Q = 3),
    each counted with the class number of the corresponding order times a
    local factor at every prime q | N/Q, namely the number of Gamma_0(q)-
    structures preserved:  1 + (d/q)  (Kronecker), except that for the
    conductor-2 order of discriminant -4Q (Q = 3 mod 4) the local factor at
    q = 2 is 2.  Validated: genus-0 quotients, Ogg's genus-0 list for
    X_0(N)^*, and the Monster-prime characterisation of g(X_0(p)^+) = 0."""
    if N % Q != 0 or gcd(Q, N // Q) != 1:
        raise ValueError(f'{Q} is not an exact divisor of {N}')
    M = N // Q
    others = list(primefact(M)) if M > 1 else []
    def loc(disc, two_factor=None):
        f = 1
        for q in others:
            if q == 2 and two_factor is not None:
                f *= two_factor
            else:
                f *= 1 + quad_rec(disc, q)
        return f
    if Q == 2:
        return _hw(-4) * loc(-4) + _hw(-8) * loc(-8)
    if Q == 3:
        return _hw(-3) * loc(-3) + _hw(-12) * loc(-12, two_factor=2)
    if Q == 4:
        return _hw(-4) * loc(-4)
    if Q % 4 == 3:
        return _hw(-4 * Q) * loc(-4 * Q, two_factor=2) + _hw(-Q) * loc(-Q)
    return _hw(-4 * Q) * loc(-4 * Q)

def genus_x0_plus(N:int)->int:
    """Genus of X_0(N)^+ = X_0(N)/w_N (Fricke quotient), N >= 2."""
    g0 = genus_x0(N)
    nu = al_fixed_points(N, N)
    g2 = 2 * (g0 - 1) - nu + 4          # 2g0 - 2 = 2(2g+ - 2) + nu
    assert g2 % 4 == 0, (N, g0, nu)
    return g2 // 4

def genus_x0_star(N:int)->int:
    """Genus of X_0(N)^* = X_0(N)/<all w_Q>, N squarefree with r prime factors:
    Riemann-Hurwitz for the (Z/2)^r quotient,
        2g0 - 2 = 2^r (2g* - 2) + sum_{Q | N, Q > 1} nu(w_Q)."""
    pf = primefact(N)
    if any(e > 1 for e in pf.values()):
        raise NotImplementedError('genus_x0_star: only squarefree N')
    r = len(pf)
    g0 = genus_x0(N)
    nu = sum(al_fixed_points(Q, N) for Q in divisors(N) if Q > 1)
    num = 2 * g0 - 2 - nu + 2 ** (r + 1)     # = 2^r (2 g*)
    assert num % (2 ** (r + 1)) == 0, (N, g0, nu)
    return num // (2 ** (r + 1))


def _selftest_genus(verbose=True):
    """Sanity battery for the genus functions.  Returns True if all pass."""
    from nt import primesBetween
    ok = True
    # genus of X_0(N) for small N (standard table)
    known = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:1,12:0,13:0,14:1,15:1,
             16:0,17:1,18:0,19:1,20:1,21:1,22:2,23:2,24:1,25:0,26:2,27:1,28:2,
             29:2,30:3,31:2,32:1,33:3,34:3,35:3,36:1,37:2,38:4,39:3,40:3,41:3,
             42:5,43:3,44:4,45:3,46:5,47:4,48:3,49:1,50:2,71:6,97:7,127:10}
    ok &= all(genus_x0(N) == g for N, g in known.items())
    # Ogg: g(X_0(p)^+) = 0 iff p | #Monster
    monster = [2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]
    ok &= [p for p in primesBetween(2, 200) if genus_x0_plus(p) == 0] == monster
    # Ogg / Conway-Norton: squarefree N with g(X_0(N)^*) = 0
    ogg = [2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,38,39,
           41,42,46,47,51,55,59,62,66,69,70,71,78,87,94,95,105,110,119]
    sqf = [N for N in range(2, 300) if all(e == 1 for e in primefact(N).values())]
    ok &= [N for N in sqf if genus_x0_star(N) == 0] == ogg
    # fixed-point parity: nu(w_Q) = 2 g0 + 2 mod 4
    ok &= all((al_fixed_points(Q, N) - 2 * genus_x0(N) - 2) % 4 == 0
              for N in sqf for Q in divisors(N) if Q > 1)
    if verbose:
        print('genus self-test:', 'OK' if ok else 'FAILED')
    return ok
