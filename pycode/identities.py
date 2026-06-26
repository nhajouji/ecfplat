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