
def gcd(a:int,b:int)->int:
    a = abs(a)
    b = abs(b)
    if min(a,b)==0:
        return max(a,b)
    while b % a != 0:
        r = b % a 
        b = a
        a = r
    return a

def lcm(a:int,b:int):
    d = gcd(a,b)
    if d!= 0:
        return (a*b)//d
    else:
        return a*b

def gcd_list(l):
    l = list({c for c in l if c != 0})
    if len(l)==0:
        return 0
    g = l[0]
    for x in l[1:]:
        g = gcd(g,x)
        if g == 1:
            return 1
    return g

def axby(ab:tuple[int])->tuple[int]:
    a, b = ab
    x0, y0, x1, y1 = 0, 1, 1, 0
    while a != 0:
        r = b % a
        q = b // a
        x2 = q*x1 + x0
        x0 = x1
        x1 = x2
        y2 = q*y1 + y0
        y0 = y1
        y1 = y2
        b = a
        a = r
    return x0,y0

def hall_multiplier(l:int,m:int)->int:
    if min(l,m)< 1:
        return 'Not defined'
    l0 = gcd(l,m)
    l1 = m//l0
    multiplier = 1
    while gcd(l0,l1)>1:
        g = gcd(l0,l1)
        l0*= g
        l1 = l1//g
        multiplier*=g
    return multiplier

# Mod p

def mod_sfd(a,m):
    a = a % m
    if 2*a > m:
        a-=m
    return a

## Prime factorization

def primefact(n:int)->dict:
    if n == 0:
        return {0:1}
    pf = {}
    if n < 0:
        n = abs(n)
        pf[-1]=1
    for p0 in [2,3]:
        if n % p0 == 0:
            pf[p0] = 0
            while n % p0 == 0:
                n = n//p0
                pf[p0]+=1
    e = -1
    p = 5
    while p*p <= n:
        if n % p == 0:
            pf[p] = 0
            while n % p == 0:
                n = n//p
                pf[p]+=1
        p+=3+e
        e*=-1
    if n>1:
        pf[n]=1
    return pf

def find_prim_root(p):
    mg_ord = p-1
    mg_max_facs = [mg_ord//m for m in primefact(mg_ord)]
    for a in range(2,p):
        aks = [pow(a,m,p) for m in mg_max_facs]
        if min(aks)>1:
            return a
    return 1

def primeQ(n:int)->bool:
    pfn = primefact(n)
    return len(pfn)==1 and max(pfn.values())==1

def pf_to_int(pf:dict)->int:
    n = 1
    for p in pf:
        n*=(p**pf[p])
    return n

def pf_to_divisors(pf:dict)->list:
    divs = [1]
    for p in pf:
        divs = [d*(p**e) for d in divs for e in range(pf[p]+1)]
    divs.sort()
    return divs

def divisors(n):
    return pf_to_divisors(primefact(n))

def no_odd_prime_facs(d):
    return len([p for p in primefact(abs(d)) if p % 2 == 1])

def quad_gcd(a1:int,a2:int)->int:
    pf2 = primefact(a2)
    a2rt = pf_to_int({p:pf2[p]//2 for p in pf2})
    return gcd(a1,a2rt)


def ap_to_lm(ap:tuple[int])->tuple[int]:
    a, p = ap
    #min poly of frob is x^2 - ax + p
    #trace of frob - 1 is a-2 and norm is p+a+1
    t = a - 2
    n = p - a + 1
    m = quad_gcd(t,n)
    l = n//(m**2)
    if n!= l*(m**2):
        return 'Something went wrong'
    return (l,m)

########################
# Quadratic characters #
########################
def quad_rec(a,p):
    if p == 2:
        if a % 2 == 0:
            return 0
        elif a % 8 == 1 or a % 8 == 7:
            return 1
        else:
            return - 1
    pwr = pow(a,p//2,p)
    return pwr - ((pwr+1)//p)*p

def find_nonsquare(p):
    """The largest non-residue in [1, p-1] (used to twist a model to a chosen signature)."""
    d = p - 1
    while quad_rec(d, p) == 1:
        d -= 1
    return d

def jacobi_symbol(d:int,n:int):
    pfn = primefact(n)
    s = 1
    for p in pfn:
        s*=(quad_rec(d,p)**pfn[p])
    return s

def gen_quad_symb(a,b):
    pf_b = primefact(b)
    if len(pf_b) == 1 and max(pf_b.values()) == 1:
        return quad_rec(a,b)
    s = 1
    for p in pf_b:
        if pf_b[p] % 2 == 1:
            s *= quad_rec(a,p)
    return s
# Applications

def class_no_formula(d):
    # We pull out the conductor, if there is 1, and assume d is a fundamental discriminant
    d, c = discfac(d)
    # We compute class number of d
    md = 1
    if d >-5:
        if d == -3:
            md = 3
        elif d == -4:
            md = 2
    qsum = sum([jacobi_symbol(d,n) for n in range(1,(abs(d)+1)//2)])
    den =2-jacobi_symbol(d,2)
    return (md*qsum*twisted_phi(d,c))//den


## Generating prime lists for testing

def primesBetween(a:int,b:int)->list[int]:
    a = max(2,a)
    if b < a:
        return []
    elif b == 2:
        return [2]
    elif b == 3:
        if a == 2:
            return [2,3]
        else:
            return [3]
    m = b
    # The following code is going to compute the set of primes < m
    # First, we create a list that will keep track of the primes-
    # the set contains m+1 elements, and for integers k <= m,
    # cands[k] = 0 records the fact that k is composite.
    # To begin, every even element other than 2 is set to 0,
    # to record the fact that 2 is the only even prime.
    cands = [0,0,1]+[(i % 2) for i in range(3,m+1)]
    # Next, we record that multiples of 3 other than 3 are composite.
    # Note that 6 is already known to be composite, 
    # so we can start at 9.
    for i in range(9,m+1,3):
        cands[i] = 0
    # Now we start going through the list of candidates to find primes p,
    # and updating the list by recording that all multiples of the prime,
    # which lie between p^2 and m inclusive, are composite.
    # We can stop once p^2 > m, since all multiples of p^2 will exceed m.
    # We start at p = 5, and we will only check odd numbers
    # that are not multiples of 3.
    # To avoid the multiples of 3, we will jump ahead by either 2 or 4
    # at each step. The number e records whether we need to jump by 2 or 4
    # at each step.
    e = -1
    p = 5
    while p**2 <= m:
        # This assumes p is prime. This is the case when p = 5,
        # and will be the case at the end of the loop when p is redefined.
        # As mentioned above, we start by recording that all multiples of
        # p between p^2 and m are composite.
        for pm in range(p**2,m+1,p):
            cands[pm]=0
        # Now we look for the next prime.
        # We alternate adding 2 and 4 to p and then checking whether
        # p is prime by checking if cands[p] is 0 or not.
        p+=3+e
        e*=-1
        while cands[p] == 0:
            p+=3+e
            e*=-1
    # When the loop ends, cands[p] = 0 if and only if p is 0,
    # so we can obtain the set of primes:
    primes = [p for p in range(a,m+1) if cands[p] == 1]
    return primes

def valuation(a:int,p:int):
    n = 0
    while a % p == 0:
        a = a//p
        n+=1
    return n

def discfac(d):
    if d == 0:
        return (0,0)
    elif d % 4 > 1:
        return (d,1)
    m = 1
    s = d//abs(d)
    d*=s
    while d % 4 == 0:
        m*=2
        d = d//4
    while d % 9 == 0:
        m*=3
        d = d//9
    r = 5
    e = -1
    while r*r <= d:
        r2 = r*r
        while d % r2 == 0:
            d = d//r2
            m*=r
        r+=3+e
        e*=-1
    d*=s
    if d % 4 > 1:
        return d*4, m//2
    else:
        return d, m 
    


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

## CRT

def crt_pair(am1:tuple[int],am2:tuple[int])->tuple[int]:
    a1,m1 = am1
    a2,m2 = am2
    if gcd(m1,m2)>1:
        return 'Check moduli'
    m12 = m1*m2
    a12 = (a1 + m1*((a2-a1)*pow(m1,-1,m2))) % m12
    if 2*a12 > m12:
        a12-=m12
    return (a12,m12)

def crt_list(amlist:list[tuple[int,int]]):
    if len(amlist)== 0:
        raise ValueError('Need at least one pair')
    am0 = amlist[0]
    for am1 in amlist[1:]:
        am0 = crt_pair(am0,am1)
    return am0

## Root of unity
def get_rou_mod(m,p):
    m0 = gcd(m,p-1)
    k = (p-1)//m0
    for x in range(2,p-1):
        xk = pow(x,k,p)
        xkg = [pow(xk,i,p) for i in range(m0)]
        if len(set(xkg))== m0:
            return xk
    return 0


def int_sqrt(n:int)->int:
    if not isinstance(n,int) or n<0:
        raise ValueError('n must be a nonnegative integer')
    if n < 2:
        return n
    r = 1
    while r**2 + 1 < n:
        b = 2
        while (r+b)**2 < n:
            b*=2
        r+= (b//2)
    return r


#####################################
# Isogeny kernel extension degrees  #
#####################################

def mult_order_mod(x:int, m:int)->int:
    """Multiplicative order of x in (Z/m)^* (x must be a unit mod m)."""
    x %= m
    if gcd(x, m) != 1:
        raise ValueError(f'{x} is not a unit mod {m}')
    k, cur = 1, x
    while cur != 1:
        cur = (cur * x) % m
        k += 1
    return k


def sqrt_mod(d:int, l:int):
    """A square root of d mod l (l prime), or None if d is a non-residue.
    Brute force -- intended for small l."""
    d %= l
    for r in range(l):
        if (r * r) % l == d:
            return r
    return None


def _order_in_fl2(a:int, p:int, l:int)->int:
    """Order of t (a Frobenius eigenvalue) in F_l[t]/(t^2 - a t + p), the inert
    case.  t^2 = a t - p, so multiplication reduces explicitly -- pure mod-l
    arithmetic, no field object needed.  Result divides l^2 - 1."""
    def mul(u, v):
        c0 = (u[0]*v[0]) % l
        c1 = (u[0]*v[1] + u[1]*v[0]) % l
        c2 = (u[1]*v[1]) % l                 # coefficient of t^2 = a t - p
        return ((c0 - c2*p) % l, (c1 + c2*a) % l)
    one, t = (1, 0), (0, 1)
    k, cur = 1, t
    while cur != one:
        cur = mul(cur, t)
        k += 1
    return k


def frob_ext_degrees(a:int, p:int, l:int)->dict:
    """Extension degrees needed to compute the l-isogenies of an ordinary curve.

    For an ordinary elliptic curve of trace a over F_p (l prime, l != p), Frobenius
    pi acts on E[l] = (Z/l)^2 with characteristic polynomial x^2 - a x + p (mod l).
    An l-isogeny kernel is a pi-eigenline; a generator P of the eigenline for
    eigenvalue lambda satisfies pi(P) = lambda*P, so P is defined over F_{p^k} with
        k = ord(lambda)   (multiplicative order of the eigenvalue).

    Returns {'kind', 'eigenvalues', 'degrees', 'min_degree'} where kind is:
        'split'    - two distinct eigenvalues in F_l (l splits: two horizontal
                     l-isogenies, the +-x_l directions),
        'ramified' - a repeated eigenvalue,
        'inert'    - eigenvalues in F_{l^2} (no horizontal l-isogeny); the single
                     degree is the order in F_{l^2}^*, dividing l^2 - 1.
    For split/ramified, eigenvalues and degrees are aligned lists over F_l, and
    min_degree is the cheapest direction to compute.
    """
    if p % l == 0:
        raise ValueError(f'l = {l} must differ from the characteristic p = {p}')
    a0, p0 = a % l, p % l
    disc = (a0 * a0 - 4 * p0) % l
    if l == 2:
        eigs = sorted({x for x in range(2) if (x*x - a0*x + p0) % 2 == 0})
        kind = 'split' if len(eigs) == 2 else ('ramified' if len(eigs) == 1 else 'inert')
    else:
        r = sqrt_mod(disc, l)
        if r is None:
            kind, eigs = 'inert', []
        else:
            inv2 = pow(2, l - 2, l)
            eigs = sorted({((a0 + r) * inv2) % l, ((a0 - r) * inv2) % l})
            kind = 'ramified' if r == 0 else 'split'
    if kind == 'inert':
        k = _order_in_fl2(a0, p0, l)
        return {'kind': 'inert', 'eigenvalues': None, 'degrees': [k], 'min_degree': k}
    degrees = [mult_order_mod(lam, l) for lam in eigs]
    return {'kind': kind, 'eigenvalues': eigs, 'degrees': degrees, 'min_degree': min(degrees)}

# ===========================================================================
# Gaussian integers Z[i] and the sum-of-two-squares algorithm
# ---------------------------------------------------------------------------
# A Gaussian integer a + b*i is stored as a tuple (a, b).  Everything here is
# exact integer arithmetic (Python bigints), so it is safe for cryptographic
# primes -- no floats anywhere, including the round-to-nearest division step.
# ===========================================================================

def _nint(a: int, n: int) -> int:
    """Nearest integer to a/n for n > 0 (ties round toward +infinity).
    Pure integer arithmetic -- no float rounding, safe for huge a."""
    return (2 * a + n) // (2 * n)


def gi_norm(z: tuple) -> int:
    a, b = z
    return a * a + b * b


def gi_conj(z: tuple) -> tuple:
    a, b = z
    return (a, -b)


def gi_add(z: tuple, w: tuple) -> tuple:
    return (z[0] + w[0], z[1] + w[1])


def gi_sub(z: tuple, w: tuple) -> tuple:
    return (z[0] - w[0], z[1] - w[1])


def gi_mul(z: tuple, w: tuple) -> tuple:
    a, b = z
    c, d = w
    return (a * c - b * d, a * d + b * c)


def gi_rounddiv(z: tuple, w: tuple) -> tuple:
    """Closest Gaussian integer to z / w  (w != 0).

    This is the geometry-of-numbers step: writing z in the orthogonal basis
    (w, i*w) of the lattice Z[i]*w, the closest multiple of w to z is obtained
    by rounding each coordinate <z, w>/N(w) and <z, i*w>/N(w) to the nearest
    integer.  Concretely both come from z * conj(w) / N(w)."""
    n = gi_norm(w)
    a, b = gi_mul(z, gi_conj(w))          # z * conj(w) = <z,w> + <z,i w> i
    return (_nint(a, n), _nint(b, n))


def gi_divmod(z: tuple, w: tuple) -> tuple:
    """Division with remainder in Z[i]: z = q*w + r with N(r) <= N(w)/2 < N(w)."""
    q = gi_rounddiv(z, w)
    r = gi_sub(z, gi_mul(q, w))
    return q, r


def gi_gcd(z: tuple, w: tuple, path: list = None) -> tuple:
    """Euclidean algorithm in Z[i].  If `path` is a list, each remainder is
    appended to it (useful for step counts and for drawing the descent)."""
    while gi_norm(w) != 0:
        z, w = w, gi_divmod(z, w)[1]
        if path is not None:
            path.append(w)
    return z


def sqrt_mod_prime(a: int, p: int):
    """A square root of a mod p (p an odd prime), or None if a is a
    non-residue.  Tonelli-Shanks -- works for cryptographic-size p."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, (t * b * b) % p, (r * b) % p
    return r


def sqrt_neg1(p: int):
    """Smallest positive r with r^2 = -1 (mod p).  Exists iff p == 2 or
    p % 4 == 1; returns None otherwise (or for p == 2, returns 1)."""
    if p == 2:
        return 1
    r = sqrt_mod_prime(p - 1, p)
    if r is None:
        return None
    return min(r, p - r)


def sos(p: int) -> tuple:
    """Write a prime p = a^2 + b^2 (requires p == 2 or p % 4 == 1), returning
    (a, b) with 0 <= a <= b.

    Method: seed z1 = 1 + r_p*i (which has norm r_p^2 + 1, a multiple of p),
    then g = gcd_{Z[i]}(z1, p).  Because Z[i] is Euclidean, (z1, p) = the prime
    ideal above p, so N(g) = p exactly.  Always succeeds -- no search."""
    if p == 2:
        return (1, 1)
    r = sqrt_neg1(p)
    if r is None:
        raise ValueError(f"{p} is not 2 and not 1 mod 4; no sum of two squares")
    g = gi_gcd((1, r), (p, 0))
    a, b = abs(g[0]), abs(g[1])
    return (min(a, b), max(a, b))


def sos_cf(p: int, path: list = None) -> tuple:
    """Sum of two squares by the Serret / Hermite continued-fraction descent
    ("trick 1"): run the ordinary integer Euclidean algorithm on (p, r_p) and
    read off the first pair of consecutive remainders that drop below sqrt(p).

    This is the small-integer, fully elementary descent -- no Gaussian
    arithmetic.  If `path` is a list, the remainder sequence is recorded."""
    if p == 2:
        return (1, 1)
    r = sqrt_neg1(p)
    if r is None:
        raise ValueError(f"{p} is not 2 and not 1 mod 4; no sum of two squares")
    import math
    a, b = p, r
    if path is not None:
        path.extend([a, b])
    root = math.isqrt(p)
    while b > root:
        a, b = b, a % b
        if path is not None:
            path.append(b)
    # b is the first remainder <= sqrt(p); (b, a % b) are the straddling pair
    c = a % b
    return (min(b, c), max(b, c))


# ===========================================================================
# Eisenstein integers Z[omega] and the p = x^2 + 3 y^2 algorithm
# ---------------------------------------------------------------------------
# omega = (-1 + sqrt(-3))/2, a primitive cube root of unity, omega^2 = -1 - omega.
# An Eisenstein integer a + b*omega is stored as (a, b).  Norm form is the
# hexagonal-lattice norm  N(a + b*omega) = a^2 - a b + b^2.
#
# Z[omega] is the ring of integers of Q(sqrt(-3)) (disc -3, class number 1) and
# is norm-Euclidean: coordinate-rounding of z*conj(w)/N(w) in the {1, omega}
# basis gives a remainder of norm <= (3/4) N(w), because the norm form is
# bounded by 3/4 on the box [-1/2,1/2]^2.  So the same gcd machinery as Z[i]
# applies, and gcd(seed, p) always has norm exactly p.
# ===========================================================================

def e_norm(z: tuple) -> int:
    a, b = z
    return a * a - a * b + b * b


def e_conj(z: tuple) -> tuple:
    """Complex conjugate: bar(a + b*omega) = a + b*omega^2 = (a-b) - b*omega."""
    a, b = z
    return (a - b, -b)


def e_add(z: tuple, w: tuple) -> tuple:
    return (z[0] + w[0], z[1] + w[1])


def e_sub(z: tuple, w: tuple) -> tuple:
    return (z[0] - w[0], z[1] - w[1])


def e_mul(z: tuple, w: tuple) -> tuple:
    """(a+b w)(c+d w) = (ac - bd) + (ad + bc - bd) w,  using w^2 = -1 - w."""
    a, b = z
    c, d = w
    return (a * c - b * d, a * d + b * c - b * d)


def e_mul_omega(z: tuple) -> tuple:
    """Multiply by omega (rotation by 120 deg): omega*(a + b w) = -b + (a-b) w."""
    a, b = z
    return (-b, a - b)


def e_rounddiv(z: tuple, w: tuple) -> tuple:
    """Closest Eisenstein integer to z / w (w != 0), by rounding the {1,omega}
    coordinates of z*conj(w)/N(w).  Valid Euclidean step (norm form <= 3/4 on
    the fundamental box)."""
    n = e_norm(w)
    a, b = e_mul(z, e_conj(w))
    return (_nint(a, n), _nint(b, n))


def e_divmod(z: tuple, w: tuple) -> tuple:
    q = e_rounddiv(z, w)
    r = e_sub(z, e_mul(q, w))
    return q, r


def e_gcd(z: tuple, w: tuple, path: list = None) -> tuple:
    """Euclidean algorithm in Z[omega]; records remainders in `path` if given."""
    while e_norm(w) != 0:
        z, w = w, e_divmod(z, w)[1]
        if path is not None:
            path.append(w)
    return z


def eisenstein_units() -> list:
    """The six units of Z[omega]: 1, omega, omega^2, and their negatives."""
    us = []
    z = (1, 0)
    for _ in range(3):
        us.append(z)
        us.append((-z[0], -z[1]))
        z = e_mul_omega(z)
    return us


def _seed_omega(p: int):
    """An Eisenstein integer of norm divisible by p (p prime, p % 3 == 1):
    solve b^2 - b + 1 = 0 mod p via b = (1 + sqrt(-3))/2, giving z = 1 + b*omega
    with N(z) = 1 - b + b^2 = 0 mod p."""
    s = sqrt_mod_prime((-3) % p, p)
    if s is None:
        return None
    inv2 = pow(2, p - 2, p)
    b = ((1 + s) * inv2) % p
    return (1, b)


def esos(p: int) -> tuple:
    """Write a prime p = a^2 - a b + b^2 (Eisenstein norm), requiring p == 3 or
    p % 3 == 1.  Returns an Eisenstein integer (a, b) of norm exactly p, via the
    Z[omega] gcd of the seed with p.  Always succeeds (Z[omega] Euclidean)."""
    if p == 3:
        return (1, 2)                      # N(1 + 2w) = 1 - 2 + 4 = 3
    seed = _seed_omega(p)
    if seed is None:
        raise ValueError(f"{p} is not 3 and not 1 mod 3; not of the form a^2-ab+b^2")
    return e_gcd(seed, (p, 0))


def x2_3y2(p: int) -> tuple:
    """Write a prime p = x^2 + 3 y^2 (requires p == 3 or p % 3 == 1), with
    x, y >= 0.  Obtain the Eisenstein rep p = c^2 - c d + d^2, then pick the unit
    associate whose second coordinate is even (always possible: one of c, d, c-d
    is even) and complete the square: p = (c - d/2)^2 + 3 (d/2)^2."""
    a, b = esos(p)
    for u in eisenstein_units():
        c, d = e_mul((a, b), u)
        if d % 2 == 0:
            y = d // 2
            x = c - y
            return (abs(x), abs(y))
    raise RuntimeError("no associate with even second coordinate (impossible)")


def frob_traces_j0(p: int) -> list:
    """The six candidate Frobenius traces for a j = 0 curve over F_p (p % 3 == 1,
    ordinary): trace of the unit-associate u*pi is 2c - d for associate (c, d).
    Returns the six traces sorted.  The actual trace of y^2 = x^3 + k is one of
    these, selected by the sextic residue class of k."""
    a, b = esos(p)
    traces = set()
    for u in eisenstein_units():
        c, d = e_mul((a, b), u)
        traces.add(2 * c - d)
    return sorted(traces)
