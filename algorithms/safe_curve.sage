# General stuff needed by multiple scripts goes here

def is_strong(r0, k0, p, N):
    """Returns r != 0 if the group order N fits the security
    requirements of factoring in N = kr with k <= k0, r >= r0
    and r prime. Returns 0 otherwise."""
    r, k = 0, 0
    # Compute r
    for i in range(1, k0+1):
        if N % i == 0 and N/i >= r0 and N/i in Primes():
            r = N/i
            k = i
    # Check that the curve isn't anomalous
    if r == 0 or r == p:
        return 0
    # Check size requirements
    if k > k0 or r < r0:
        return 0
    # Check that the multiplicative order of p in F_r is at least 20
    p_r = Mod(p, r)
    if any(p_r^i == 1 for i in range(1, 20+1)):
        return 0
    return r

def get_generator(E, k):
    """Returns a point on E of order #E(GF(p))/k."""
    while True:
        G = k*E.random_point()
        if G != 0*G:
            return G
