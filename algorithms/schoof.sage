#!/usr/bin/sage

# Implements Schoof's point counting algorithm on elliptic curves

load('safe_curve.sage')

# Curve operations
def double(P, a):
    if P == None:
        return None
    x = P[0]
    y = P[1]
    if y == 0:
        return None
    slope = (3*x*x + a)*(2*y)**(-1)
    x_res = slope**2 - 2*x
    y_res = slope*(x - x_res) - y
    return (x_res, y_res)

def add(P, Q, a):
    if P == None:
        return Q
    if Q == None:
        return P
    if P == Q:
        return double(P, a)
    x_p, y_p = P[0], P[1]
    x_q, y_q = Q[0], Q[1]
    if x_p == x_q:
        return None
    slope = (y_p - y_q)/(x_p - x_q)
    v = (y_p*x_q - y_q*x_p)/(x_q - x_p)
    x_res = slope*slope - x_p - x_q
    y_res = -slope*x_res - v
    return (x_res, y_res)

def mul(k, P, a):
    if k == 0:
        return None
    if k % 2 == 0:
        return mul(k//2, double(P, a), a)
    return add(P, mul(k//2, double(P, a), a), a)

def frob(P, p):
    """Implements the p-th Frobenius map."""
    return (P[0]^p, P[1]^p)

def get_quotient_ring(a, b, p, l):
    """Returns F_q[x, y]/(psi_l, y^2 - f) where psi_l is the l-th
    division polynomial of E, and f such that y^2 = f(x)."""
    R.<xbar, ybar> = PolynomialRing(GF(p))
    E = EllipticCurve(GF(p), [a, b])
    psi_l = E.division_polynomial(l, xbar)
    f = xbar^3 + a*xbar + R(b)
    return R.quotient([psi_l, ybar^2 - f], ['x', 'y'])

def schoof(a, b, p):
    """Returns #E(F_p) where E is the elliptic curve over F_p defined by
    y^2 = x^3 + ax + b, using Schoof's algorithm."""
    A = 2
    l = 3
    primes = [2]
    n_mods = []

    # Special case for l=2
    w = PolynomialRing(GF(p), 'w').gen()
    f = w^3 + a*w + b
    if f.gcd(w^p - w) == 1:
        n_mods = [1]
    else:
        n_mods = [0]

    while A < 4*sqrt(p):
        R = get_quotient_ring(a, b, p, l)
        x, y = R.gens()

        frob_p = frob((x, y), p)
        frob_p2 = frob(frob_p, p)

        # (x^p2, y^p2) + [p % l](x, y)
        tmp = mul(p % l, (x, y), a)
        lhs = add(frob_p2, tmp, a)

        for n in range(l):
            if n > 0:
                rhs = add(rhs, frob_p, a)
            else:
                rhs = None

            if lhs == rhs:
                n_mods.append(n)
                primes.append(l)
                A *= l
                break
        l = next_prime(l)

    trace = crt(n_mods, primes)

    if trace > 2*sqrt(p):
        trace -= A

    return p + 1 - trace

def division_polynomial(n, p, a4, a6):
    """Returns the n-th division polynomial of the curve
    E/F_p: y^2 = x^3 + a4*x + a6."""
    R.<x, y> = GF(p)[]
    psis = [R(0), R(1), 2*y, 3*x^4 + 6*a4*x^2 + 12*a6*x - a4^2,
        4*y*(x^6 + 5*a4*x^4 + 20*a6*x^3 - 5*a4^2*x^2 - 4*a4*a6*x\
                -8*a6^2 - a4^3)]
    def psi(n):
        if n <= 4:
            return psis[n]
        m = n // 2
        if n % 2 == 0:
            return psi(m)*(psi(m+2)*psi(m-1)^2-psi(m-2)*psi(m+1)^2)/(2*y)
        return psi(m+2)*psi(m)^3-psi(m-1)*psi(m+1)^3
    return psi(n)
