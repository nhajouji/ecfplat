from nt import get_rou_mod, gcd, gcd_list, hall_multiplier, axby

############
# Matrices #
############

def shape(m:list[list[int]]):
    d2s = list({len(r) for r in m})
    if len(d2s)>1:
        raise ValueError('Not a matrix')
    d2 = d2s[0]
    return (len(m),d2)

def dot(v1:list,v2:list):
    if len(v1)!=len(v2):
        raise ValueError('Vectors must be of equal dimension')
    return sum({x*y for x,y in zip(v1,v2)})


def mat_mul(m1,m2):
    d0,d1 = shape(m1)
    d2,d3 = shape(m2)
    if d1 != d2:
        raise ValueError('Incompatible matrices')
    mprod = [[_ for _ in range(d3)] for _ in range(d0)]
    for i, r in enumerate(m1):
        for j in range(d3):
            cj = [r2[j] for r2 in m2]
            mprod[i][j]=dot(r,cj)
    return mprod


def check_square_matrix(m)->bool:
    s= {len(r) for r in m}
    if len(s)>1:
        return False
    l0 = min(s)
    if l0!= len(m):
        return False
    return True

def dotvec(v1,v2):
    if len(v1)!=len(v2):
        raise ValueError('Lengths must be equal')
    s = 0
    for i,x in enumerate(v1):
        s+=x*v2[i]
    return s

class IntegerSquareMatrix:
    def __init__(self,data:list[list[int]]):
        if type(data)!= list or type(data[0])!= list:
            raise ValueError('Input must be list of lists')
        if len(data)!=len(data[0]):
            raise ValueError('Input must be square matrix')
        self.mat = data
        self.dim = len(data)

    def __str__(self):
        m = self.mat
        s = '['
        for r in m[:-1]:
            s+=str(r)+',\n '
        s+=str(m[-1])+']'
        return s
    
    def __repr__(self):
        return str(self)
    
    def __eq__(self, other):
        if not isinstance(other, IntegerSquareMatrix):
            return False
        return self.mat == other.mat
                
    def mvec(self,v:list):
        if len(v)!=self.dim:
            raise ValueError('Wrong dimension')
        return [dot(r,v) for r in self.mat]
    
    def gcd(self):
        g_rows = [gcd_list(r) for r in self.mat]
        return gcd_list(g_rows)
    
    def gcdfac(self):
        m = self.mat
        d = self.gcd()
        if d < 2:
            return self,d
        return IntegerSquareMatrix([[a//self.gcd() for a in r] for r in self.mat]),self.gcd()
    
    
    def transpose(self):
        mat = self.mat
        return IntegerSquareMatrix([[r[i] for r in mat] for i in range(len(mat))])
    
    def trace(self):
        mat = self.mat
        d = self.dim
        return sum([mat[i][i] for i in range(d)])
    
    def det(self):
        if self.d!= 2:
            raise ValueError('Not yet implemented')
        mat = self.mat
        a,b,c,d = mat[0][0],mat[0][1],mat[1][0],mat[1][1]
        return a*d-b*c
        
    
    def __rmul__(self,n:int):
        if not isinstance(n,int):
            return self*n
        mat = self.mat
        return IntegerSquareMatrix([[x*n for x in r] for r in mat])
    
    def __neg__(self):
        return (-1)*self
    
    def __add__(self,other:'IntegerSquareMatrix'):
        if isinstance(other,int):
            one = IntegerSquareMatrix([[int(i==j) for i in range(self.dim)]
                                       for j in range(self.dim)])
            other = other*one
        if not isinstance(other,IntegerSquareMatrix):
            raise ValueError('Can only add matrices')
        m = self.dim
        if other.dim != m:
            raise ValueError('Dimensions must be equal')
        mat1 = self.mat
        mat2 = other.mat
        matsum = [[0 for i in range(m)] for j in range(m)]
        for i, r in enumerate(mat1):
            for j, x in enumerate(r):
                matsum[i][j]=x+mat2[i][j]
        return IntegerSquareMatrix(matsum)
    
    def __mul__(self,other:'IntegerSquareMatrix'):
        if isinstance(other,int):
            one = IntegerSquareMatrix([[int(i==j) for i in range(self.dim)]
                                       for j in range(self.dim)])
            other = other*one
        if not isinstance(other,IntegerSquareMatrix):
            raise ValueError('Can only multiply square matrices')
        m = self.dim
        if other.dim != m:
            raise ValueError('Dimensions must be equal')
        mat1 = self.mat
        mat2 = (other.transpose()).mat
        matprod = [[0 for i in range(m)] for j in range(m)]
        for i, r1 in enumerate(mat1):
            for j, r2 in enumerate(mat2):
                matprod[i][j] = dotvec(r1,r2)
        return IntegerSquareMatrix(matprod)
    
    def __pow__(self,n:int):
        if not isinstance(n,int) or n <0:
            raise ValueError('Can only raise to nonnegative integer powers')
        m = self.dim
        mn = IntegerSquareMatrix([[int(i==j) for i in range(m)] for j in range(m)])
        m2k = self
        while n > 0:
            if n % 2 == 1:
                mn = mn*m2k
            m2k = m2k * m2k
            n = n //2
        return mn
    
    
    
    def __sub__(self,other:'IntegerSquareMatrix')->'IntegerSquareMatrix':
        return self+(-1)*other
    
    def polyeval(self,coefs_x0_to_xn:list[int]):
        one = self**0
        polymat = 0*one
        for c in coefs_x0_to_xn[::-1]:
            polymat = self*polymat+c*one
        return polymat



def poly2mat(coefs_x0to_xn1:list[int]):
    dim = len(coefs_x0to_xn1)
    mat = [[int(i==j+1) for j in range(dim)] for i in range(dim)]
    for i in range(dim):
        mat[i][-1]-=coefs_x0to_xn1[i]
    return IntegerSquareMatrix(mat)

def monicpolyremainder(moniccoefs_x0_to_xn,qcoefs):
    mat = poly2mat(moniccoefs_x0_to_xn[:-1])
    polyatmat=mat.polyeval(qcoefs)
    return [r[0] for r in polyatmat.mat]


###

def remove_lead_zeros(l,char=0):
    if char!=0:
        l = [c% char for c in l]
    if len(l)==0:
        return [0]
    while len(l)>1 and l[-1]==0:
        l = l[:-1]
    return l

class Polynomial:
    def __init__(self,coefs_0_to_n:list[int],char=0):
        self.coefs = remove_lead_zeros(coefs_0_to_n,char)
        self.char = char
        self.lc = (self.coefs)[-1]
        self.deg = len(self.coefs)-1-int(self.lc==0)
        self.const = (self.coefs)[0]
    def __str__(self):
        coefs = self.coefs
        if len(coefs)==1:
            return str(coefs[0])
        s = ''
        for i, c in enumerate(coefs):
            if c != 0:
                if len(s)>0:
                    if s[0]=='-':
                        s = f'{c} x^{i}'+s
                    else:
                        s = f'{c} x^{i}+'+s
                else:
                    s = f'{c} x^{i}'
        if len(s)>0:
            return s
        else:
            return '0'
    def __repr__(self):
        return str(self)
    def __add__(self,other:'Polynomial'):
        if isinstance(other,int):
            other = Polynomial([other],self.char)
        if not isinstance(other,Polynomial) or self.char!=other.char:
            raise ValueError('Can only add polynomials of equal char')
        cs1 = self.coefs
        cs2 = other.coefs
        l1 = len(cs1)
        l2 = len(cs2)
        cs1+=max(l2-l1,0)*[0]
        cs2+=max(l1-l2,0)*[0]
        cssum = remove_lead_zeros([cs1[i]+cs2[i] for i in range(max(l1,l2))],self.char)
        return Polynomial(cssum,self.char)
    def __rmul__(self,n:int):
        cs1 = self.coefs
        csn = [n*c for c in cs1]
        if self.char!=0:
            csn = [c % self.char for c in csn]
        return Polynomial(csn,self.char)
    def __sub__(self,other:'Polynomial'):
        return self+(-1)*other
    def __mul__(self,other:'Polynomial'):
        if not isinstance(other,(Polynomial,int)):
            raise ValueError('Can only multiply by polynomials or integers')
        if isinstance(other,int):
            return Polynomial([c*other for c in self.coefs],self.char)
        if isinstance(other,Polynomial) and self.char != other.char:
            raise ValueError('characteristics must be equal')
        cs1 = self.coefs
        cs2 = other.coefs
        cs12 = [0 for _ in range(len(cs1)+len(cs2)-1)]
        for i, c1 in enumerate(cs1):
            for j, c2 in enumerate(cs2):
                cs12[i+j]+=c1*c2
        if self.char > 0:
            cs12 = [c% self.char for c in cs12]
        return Polynomial(cs12,self.char)
    def __pow__(self,n):
        if n < 0:
            raise ValueError('Exponent must be nonnegative')
        xn = Polynomial([1],self.char)
        x2n = self
        while n >0:
            if n % 2 == 1:
                xn = xn* x2n
            n = n//2
            x2n =x2n*x2n
        return xn
    
    def monic_associate(self):
        if self.lc <0:
            self= (-1)*self
        if self.deg == -1:
            return self
        elif self.lc == 1:
            return self
        elif self.char>0 and self.lc!=0:
            lcinv = pow(self.lc,-1,self.char)
            return lcinv *self
        else:
            g = gcd_list(self.coefs)
            return Polynomial([c//g for c in self.coefs])
        
    def __mod__(self,den:'Polynomial')->'Polynomial':
        if not isinstance(den,Polynomial) or den.char != self.char:
            raise ValueError('Can only mod polynomials of equal char')
        if den.deg <0:
            raise ZeroDivisionError('Denominator must be nonzero')
        if den.deg == 0:
            return 0*self
        den0 = den.monic_associate()
        rem0 = self.monic_associate()
        px = Polynomial([0,1],self.char)
        while rem0.deg >= den0.deg:
            lcr = rem0.lc
            lcd = den0.lc
            dd = rem0.deg-den0.deg
            rem_new = ((lcd*rem0-(px**dd)*lcr*den0)).monic_associate()
            # If something went wrong, this avoids infinite loops
            if rem_new.deg < rem0.deg:
                rem0 = rem_new
            else:
                return rem0
        return rem0
    
    def mod(self,p:int):
        if self.char>0 and self.char!=p:
            raise ValueError('Already in positive characteristic')
        return PolyFp(self.coefs,char=p)
    
    def eval(self,x):
        p = self.char
        coefs = self.coefs
        evx = 0*x
        one = x**0
        if p > 0:
            for c in coefs[::-1]:
                evx = (evx*x+c*one)%p
            return evx
        else:
            for c in coefs[::-1]:
                evx = (evx*x+c*one)
            return evx
        
        

    def dx(self):
        if len(self.coefs)==1:
            return 0*self
        cs = self.coefs
        dcs = [i*c for i,c in enumerate(cs)][1:]
        return Polynomial(dcs,char=self.char)




def poly_gcd(poly1,poly2):
    r = poly2 % poly1
    while r.deg>=0:
        poly2 = poly1
        poly1 = r
        r = poly2 % poly1
    return poly1


##########
# F_p[x] #
##########

class PolyFp(Polynomial):
    def __init__(self,coefs,char):
        super().__init__(coefs,char)

    def fp_factor(self):
        p = self.char
        xp = PolyFp([0,1],p)
        fpp = xp**p - xp
        return poly_gcd(self,fpp)
    
    def no_fp_roots(self):
        return self.fp_factor().deg
    
    
    def find_roots_BrFo(self):
        p = self.char
        return [x for x in range(p) if self.eval(x)==0]






#######
# ROI #
#######

class QuadraticInteger:
    def __init__(self,elementcoefs:tuple[int],ringpoly:tuple[int]):
        self.poly = ringpoly
        self.matrix = poly2mat(ringpoly).polyeval(elementcoefs)
        self.vec = (self.matrix.transpose().mat)[0]
    def __str__(self):
        x,y = self.vec
        c0,c1 = self.poly
        return f'({x})+({y} X) mod (X^2 + {c1} X + {c0})'
    
    def __repr__(self):
        return str(self)
    
    def __add__(self,other):
        if not isinstance(other,QuadraticInteger):
            if isinstance(other,int):
                rp = self.poly
                return self+QuadraticInteger((other,0),rp)
            else:
                raise TypeError('Incompatible types')
        if self.poly != other.poly:
            raise TypeError('Use same polynomial')
        vec1 = self.vec
        vec2 = other.vec
        vecsum = (vec1[0]+vec2[0],vec1[1]+vec2[1])
        return QuadraticInteger(vecsum,self.poly)
    def __mul__(self,other):
        if not isinstance(other,QuadraticInteger):
            if isinstance(other,int):
                rp = self.poly
                return self*QuadraticInteger((other,0),rp)
            else:
                raise TypeError('Incompatible types')
        if self.poly != other.poly:
            raise TypeError('Use same polynomial') 
        matprod = (self.matrix)*(other.matrix)
        vecprod = (matprod.transpose().mat)[0]
        return QuadraticInteger(vecprod,self.poly)
    def __rmul__(self,n:int):
        x,y = self.vec
        poly = self.poly
        return QuadraticInteger((n*x,n*y),self.poly)
    
    def __sub__(self,other):
        return self+(-1)*other
    
    def trace(self):
        return self.matrix.trace()
    
    def conjugate(self):
        t = self.trace()
        return QuadraticInteger((t,0),self.poly)-self
    
    def norm(self):
        return ((self *self.conjugate()).vec)[0]
    
    def dot(self,other):
        if not isinstance(other,QuadraticInteger) or self.poly != other.poly:
            raise TypeError('Use same polynomial')
        return (self*other.conjugate()).trace()
    

#####################
# Gaussian Integers #
#####################

def smallest_rem(a,p):
    a0 = a % p
    a1 = a0-p
    if abs(a1)<a0:
        return a1
    else:
        return a0
    
def smallest_q(a,p):
    return (a-smallest_rem(a,p))//p


def improve_approx_gauss(xy,p):
    x,y = xy
    k = p//(abs(y))
    x1 = (x*k)%p
    x2 = (x*(k+1))%p
    y1 = p-k*abs(y)
    y2 = (k+1)*abs(y)-p
    v1 = (x1,y1)
    v2 = (x2,y2)
    n1 = x1**2+y1**2
    n2 = x2**2+y2**2
    if n1 == p or n1 <= n2:
        return v1
    else:
        return v2


class GaussianInteger:
    def __init__(self,ab:tuple[int,int]):
        self.vec = ab
        self.re = ab[0]
        self.im = ab[1]
        self.qi = QuadraticInteger(ab,(1,0))
        self.tr = (self.qi).trace()
        self.nrm = (self.qi).norm()
    def __str__(self):
        if self.nrm == 0:
            return '0'
        a,b = self.vec
        if b == 0:
            return str(a)
        if a == 0:
            return str(b)+'i'
        if b > 0:
            return str(a)+'+ i'+str(b)
        else:
            return str(a)+'- i'+str(abs(b))
    def __repr__(self):
        return str(self)
    def __add__(self,other):
        x1,y1 = self.vec
        x2,y2 = other.vec
        return GaussianInteger((x1+x2,y1+y2))
    def __sub__(self,other):
        x1,y1 = self.vec
        x2,y2 = other.vec
        return GaussianInteger((x1-x2,y1-y2))
    def __rmul__(self,n:int):
        x1,y1 = self.vec
        return GaussianInteger((n*x1,n*y1))
    def __mul__(self,other):
        x1,y1 = self.vec
        x2,y2 = other.vec
        return GaussianInteger((x1*x2-y1*y2,x1*y2+x2*y1))
    def conj(self):
        x,y = self.vec
        return GaussianInteger((x,-y))
    def tracedot(self,other):
        return ((self*other.conj()).tr)//2
    def __floordiv__(self,other):
        if not isinstance(other,GaussianInteger):
            raise ValueError('Can only divide by other Gaussian integers')
        n = other.nrm
        if n == 0:
            raise ZeroDivisionError('Denominator must have nonzero norm')
        v1 = other
        v2 = GaussianInteger((0,1))*v1
        c1 = smallest_q(v1.tracedot(self),v1.tracedot(v1))
        c2 = smallest_q(v2.tracedot(self),v2.tracedot(v2))
        return GaussianInteger((c1,c2))
    def __mod__(self,other):
        q = self//other
        return self - q*other


def gcd_gauss(z1:GaussianInteger,z2:GaussianInteger)->GaussianInteger:
    if z1.nrm == 0:
        return z2
    while (z2%z1).nrm != 0:
        r = z2 % z1
        z2 = z1
        z1 = r
    return z1


def sum_of_squares_rep(p):
    if p % 4 == 3:
        return None
    if p == 2:
        return GaussianInteger((1,1))
    # r0 is a square root of -1 mod p
    # We have to choose between r0 and p-r0;
    # We pick whichever is smallest
    r0 = min(get_rou_mod(4,p),p-get_rou_mod(4,p))
    # Our first guess is 1 +r0.
    # If p = 1+n^2, then n is necessarily r0
    z1 = GaussianInteger((1,r0))
    if z1.nrm == p:
        return z1
    # We can construct a second vector using the
    # other square root of -1. This will definitely
    # not work, but it could help us improve first guess.
    z0 = GaussianInteger((1,p-r0))
    while z1.nrm !=p and z1.nrm < z0.nrm:
        if gcd(z0.nrm,z1.nrm)< z1.nrm:
            z01 = gcd_gauss(z0,z1)
            z01c = gcd_gauss(z0,z1.conj())
            if z01.nrm % p == 0:
                if z01.nrm == p:
                    return z01
                z1 = z01
            elif z01c.nrm % p == 0:
                if z01c.nrm == p:
                    return z01c
                z1 = z01c
        if z1.nrm == p:
            return z1
        z2 = GaussianInteger(improve_approx_gauss(z1.vec,p))
        z0 = z1
        z1 = z2
    return z1