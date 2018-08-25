#!/usr/bin/env sage

# Implements multiple algorithms to solve DLPs on elliptic curve groups

from math import sqrt, ceil

import matplotlib.pyplot as plt

load('rd_approach.sage')

def baby_steps(P, Q, n = 0):
    """Searches for m such that [m]P = Q using the baby-steps giant-steps
    algorithm."""
    n = P.order()
    N = Integer(ceil(sqrt(n)))
    R = -(N*P)
    b_steps = [k*P for k in range(1, N+1)]
    g_steps = [Q+k*R for k in range(1, N+1)]
    for i, elem in enumerate(b_steps):
        try:
            j = g_steps.index(elem)
            return i+1 + (j+1)*N
        except:
            continue
    return None

# Auxiliary binning function for Pollard-rho
def f(R, a, b, P, Q):
    if Integer(R[1]) % 3 == 0:
        return P+R, a+1, b
    if Integer(R[1]) % 3 == 1:
        return 2*R, 2*a, 2*b
    else:
        return Q+R, a, b+1

def pollard_rho(P, Q):
    """Returns m such that Q = [m]P using the Pollard-rho algorithm."""
    n = P.order()
    R = 0*P
    R2 = P
    while R != R2:
        if R == 0*P:
            R, a, b = P, Mod(1, n), Mod(0, n)
            R2, a2, b2 = P, Mod(1, n), Mod(0, n)
        R, a, b = f(R, a, b, P, Q)
        Rt, at, bt = f(R2, a2, b2, P, Q) # temp variables
        R2, a2, b2 = f(Rt, at, bt, P, Q)
    assert b != b2, "Error: b == b2"
    return Integer((a - a2)/(b2 - b))

def polhig_hellman(G, Q):
    """Returns m such that Q = [m]G using the Polhig-Hellman algorithm."""
    n = G.order()
    factorization = factor(n)
    if len(factorization) == 1:
        p, e = factorization[0]
        x = 0
        gamma = (p^(e-1))*G
        for k in range(e):
            h = (p^(e-1-k))*((-x)*G + Q)
            d = baby_steps(gamma, h)
            x += p^k * d
        return x
    x_list, mod_list = [], []
    for p, e in factorization:
        power = p^e
        G_i = (n//power)*G
        H_i = (n//power)*P
        x_list.append(polhig_hellman(G_i, H_i))
        mod_list.append(power)
    return crt(x_list, mod_list)
