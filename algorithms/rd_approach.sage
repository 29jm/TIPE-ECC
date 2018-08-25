#!/usr/bin/env sage

# Generates curves through a naive random approach

load('safe_curve.sage')

def get_prime(r0, k0):
    """Returns a random prime p >= r0*k0, but not by much."""
    prod = r0*k0
    prime_index = randrange(0, int(ln(prod)))
    p = next_prime(prod)
    for i in range(prime_index):
        p = p.next_probable_prime()
    if not is_prime(p):
        p = next_prime(p)
    return p

def get_parameters(p):
    """Returns (a, b) such that EllipticCurve(GF(p), [a, b]) has
    non-zero discriminant. a and b are chosen at random."""
    d = 0
    while not d:
        a = randrange(0, p)
        b = randrange(0, p)
        d = (4*a**3 + 27*b**2) % p
    return (a, b)

def get_quadratic_nonresidue(p):
    """Returns d in F_p where d isn't a square in F_p."""
    while True:
        d = randrange(1, p)
        if kronecker(d, p) == -1:
            return d

def try_get_curve(E, r0, k0, N):
    """Returns a curve (E, r, k, G) if it satisfies the given requirements,
    None otherwise."""
    p = E.base_field().order()
    r = is_strong(r0, k0, p, N)
    if r != 0:
        k = Integer(N/r)
        r = Integer(r)
        G = get_generator(E, k)
        return E, r, k, G
    return None

def gen_random(r0, k0, custom_schoof=False, twist=False):
    """Returns (E, r, k, G) such that |E(F_p)| = k*r with
    k <= k0, r >= r0 and r prime.
    G is a generator of a subgroup of E(F_p) of prime order r.
    If `twist` is True, also check the twist of every generated curve."""
    p = get_prime(r0, k0)

    while True:
        a, b = get_parameters(p)
        E = EllipticCurve(GF(p), [a, b])
        if custom_schoof:
            N = schoof(a, b, p)
        else:
            N = E.cardinality()
        C = try_get_curve(E, r0, k0, N)

        if C:
            return C

        if not twist:
            continue

        Et = E.quadratic_twist()
        Nt = 2*p + 2 - N
        Ct = try_get_curve(Et, r0, k0, Nt)

        if Ct:
            return Ct
