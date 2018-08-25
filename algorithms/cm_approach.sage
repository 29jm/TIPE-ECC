# Generates curves through the complex multiplication method

# Entirely taken from the second paragraph of section 3.1 of:
#   `Efficient ephemeral elliptic curve cryptographic keys`
# available at https://eprint.iacr.org/2015/647.pdf

load('safe_curve.sage')

def get_special_square_free(r):
    """Returns a square-free integer d with -d mod 4 = 0 or 1,
    d != 1, 3 and d < r."""
    while True:
        d = randrange(2, r)
        if d != 3 and -d % 4 in [0, 1] and is_squarefree(d):
            return d

def find_uvp(d, r):
    """returns (u, p) such that there is a v verifying u^2+dv^2 = p
    where p >= r."""
    sq_r = int(sqrt(r))
    m = 4 if d % 4 == 3 else 1
    while True:
        # u >= 1 avoids divisions by zero later
        u, v = randrange(1, sq_r), randrange(0, sq_r)
        p = u^2 + d*v^2

        if u == 1 and d % 4 == 3: # Prevents anomalous curves
            continue

        if m == 4:
            if not p % 4 == 0:
                continue
            p /= 4

        if p >= r and is_pseudoprime(p):
            return u, v, p

def order_divides(E, N):
    """Returns whether #E(GF(p)) divides N."""
    return N % E.random_point().order() == 0

def gen_cm(r0, k0, dmax=200):
    """Returns (E, r, k, G) where #E(GF(p)) == rk and E is a
    secure elliptic curve according to `Generation Method of
    Elliptic Curve`, and G is an element of E of prime order r."""
    r, N = 0, 0

    while r == N == 0:
        d = get_special_square_free(dmax)
        u, v, p = find_uvp(d, r0)
        s = 1 if d % 4 == 3 else 2

        N1 = p+1 - s*u
        N2 = p+1 + s*u

        r1 = is_strong(r0, k0, p, N1)
        if r1:
            r, N = r1, N1
        else:
            r2 = is_strong(r0, k0, p, N2)
            if r2:
                r, N = r2, N2
            else:
                continue

    H = hilbert_class_polynomial(-d)
    j = H.any_root(GF(p))
    a = -27*j/(4*(j-1728))
    k = Integer(N/r)

    E = EllipticCurve(GF(p), [a, -a])

    if order_divides(E, N):
        G = get_generator(E, k)
        return E, r, k, G
    else:
        Et = E.quadratic_twist()
        G = get_generator(Et, k)
        return Et, r, k, G
