︠54588ef7-b4c9-4f83-ba69-a090491e452fs︠
# The elliptic curve E is in Weierstrass form y^2=f(x)=x^3+Ax+B

divpoly_factor = 0    # global variable for factor of the division polynomial when ZeroDivisionError's occur

# Elements of End(E[ell]) are represented as pairs (a,b*y), with a,b in Fp[x]/(h(x)), where h is the ell-th divpoly (or a factor of it, for example, the kernel polynomial of an isogeny)
# The y is implicit, but must be accounted for when applying the group law -- using the curve equation y^2=f(x) we can replace y^2 with f(x) whenever it appears (this effectively hides all the y's)

# In many of the functions below we pass in both A and f
# where f is the image of x^3+Ax+B in Fp[x]/(h(x)) -- we need both because if deg(h)<= 3 we cannot recover A from (x^3+Ax+B) mod h(x)
    
def add(P,Q,A,f):
    """add endomorphisms P and Q in End(E[ell])"""
    global divpoly_factor
    if not P: return Q
    if not Q: return P
    a1 = P[0]; b1 = P[1]; a2=Q[0]; b2=Q[1]
    if a1 == a2:
        if b1 == b2: return dbl(P,A,f)
        else: return ()
    try:
        m = (b2-b1)/(a2-a1)
    except ZeroDivisionError:
        ### given that a2-a1 is already reduced mod h, a ZeroDivisionError means that gcd(a2-a1,h) must be a nontrivial divisor g of h
        ### raise an error so that we can restart the algorithm working in a smaller quotient ring
        divpoly_factor = a2-a1
        raise
    a3 = f*m^2 -a1 - a2
    b3 = m*(a1-a3) - b1
    return (a3,b3)
    
def dbl(P,A,f):
    """double the endomorphism P in End(E[ell]) """
    global divpoly_factor
    if not P: return P
    a1 = P[0]; b1 = P[1]
    try:
        m = (3*a1^2+A) / (2*b1*f)
    except ZeroDivisionError:
        divpoly_factor = 2*b1*f
        raise
    a3 = f*m^2 - 2*a1
    b3 = m*(a1-a3) - b1
    return (a3,b3)

def neg(P):
    """ negate the endomorphism P in End(E[ell]) """
    if not P: return P
    return (P[0],-P[1])
    
def smul (n,P,A,f):
    """ compute the scalar multiple n*P in End(E[ell]) using double and add"""
    if not n: return ()
    nbits = n.digits(2)
    i = len(nbits)-2
    Q = P
    while i >= 0:
        Q = dbl(Q,A,f)
        if nbits[i]: Q = add(P,Q,A,f)
        i -= 1
    return Q
    
def mul (P,Q):
    """ compute the product (i.e. composition of) P*Q of two endomorphisms in End(E[ell]) """
    return (P[0].lift()(Q[0]),P[1].lift()(Q[0])*Q[1])
        
def trace_mod (E, ell):
    """ compute the trace of Frobenius of E modulo ell """
    FF=E.base_ring()
    q = FF.cardinality()                        # finite field FF_q
    R.<t>=PolynomialRing(FF)
    A=E.a4(); B=E.a6()                          # E: y^2 = x^3 + Ax + B
    if ell == 2:                                # t is odd iff f is irreducible
        if (t^3+A*t+B).is_irreducible(): return 1
        else: return 0
    h = E.division_polynomial(ell,t,0).monic()
    while true:
        try:
            RR.<x> = R.quotient(ideal(h))       # RR is End(E[ell]) (or a subring thereof)
            f = x^3+A*x+B
            xq = x^q; yq = f^((q-1)/2)
            pi = (xq,yq)                        # pi is the Frobenius endomorphism
            pi2 = mul(pi,pi)                    # pi2 = pi^2
            id = (x,RR(1))                      # identity aka mult-by-1 map
            Q = smul(q%ell,id,A,f)              # Q is the mult-by-q map
            S = add(pi2,Q,A,f)                  # S = pi^2 + q = t*pi
            if not S: return 0                  # if S=0 then t=0
            if S == pi: return 1                # if S=pi then t=1
            if neg(S) == pi: return -1          # if S=-pi then t=-1
            P = pi
            for t in range(2,ell-1):
                P = add(P,pi,A,f)               # P = t*pi
                if P==S: return t               # if S=P then we have found t
            print "Error, Frob satisfies no charpoly!!"
            assert false
        except ZeroDivisionError:
            h = gcd(h,divpoly_factor.lift())    # if we hit a zero divisor, start over with new h
            print "found %d-divpoly factor of degree %d"%(ell,h.degree())

def Schoof(E):
    """ compute the trace of Frobenius of E using Schoof's algorithm """
    q=E.base_ring().cardinality()
    t = 0; M=1; ell=1;
    while M <= 4*sqrt(q):
        ell = next_prime(ell)
        start = cputime()
        tell = trace_mod(E,ell)
        print "trace %d mod %d computed in %.2f secs"%(tell,ell,cputime()-start)
        a = M*M.inverse_mod(ell); b = ell*ell.inverse_mod(M)
        M *= ell
        t = (a*tell+b*t) % M
    if t >= M/2: return t-M
    else: return t

︡f8b6fb51-5aa1-4e20-a263-69733080cb46︡{"done":true}︡
︠10ff9a64-1e3b-4555-a103-d3cf4309d280s︠
%time
FF=GF(next_prime(2^80))
E=EllipticCurve([FF(314159),FF(2781828)])
t=Schoof(E)
print t,E.trace_of_frobenius()
︡7ae4ac16-946a-4ba6-a446-482eda997609︡{"stdout":"trace 1 mod 2 computed in 0.00 secs\ntrace -1 mod 3 computed in 0.02 secs\ntrace 0 mod 5 computed in 0.03 secs\ntrace 2 mod 7 computed in 0.04 secs"}︡{"stdout":"\ntrace 7 mod 11 computed in 0.19 secs"}︡{"stdout":"\nfound 13-divpoly factor of degree 6"}︡{"stdout":"\ntrace 6 mod 13 computed in 0.28 secs\ntrace 15 mod 17 computed in 0.80 secs"}︡{"stdout":"\ntrace 1 mod 19 computed in 0.75 secs"}︡{"stdout":"\ntrace 8 mod 23 computed in 1.86 secs"}︡{"stdout":"\ntrace 22 mod 29 computed in 6.39 secs"}︡{"stdout":"\ntrace 13 mod 31 computed in 9.21 secs"}︡{"stdout":"\ntrace 17 mod 37 computed in 38.06 secs"}︡{"stdout":"\n"}︡{"stdout":"-1315484487805 -1315484487805\n"}︡{"stdout":"CPU time: 57.75 s, Wall time: 57.90 s\n"}︡{"done":true}︡
︠f5a83690-c8ad-4f96-96dd-0d85472522e9i︠
︡04956c6c-0663-49e6-9e0a-d4261050c2ca︡
︠5de2e567-d14b-4c6e-af70-2cc88a5398b4︠









