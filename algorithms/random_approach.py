#!/bin/python

# P.S: This file has been thoroughly superseeded by rd_approach.sage

from random import randrange
from math import log
from hashlib import sha1
from math import sqrt
import subprocess

debug = True

def miller_rabin(n, k=10):
	if n == 2:
		return True
	if not n & 1:
		return False

	def check(a, s, d, n):
		x = pow(a, d, n)
		if x == 1:
			return True
		for i in range(s - 1):
			if x == n - 1:
				return True
			x = pow(x, 2, n)
		return x == n - 1

	s = 0
	d = n - 1

	while d % 2 == 0:
		d >>= 1
		s += 1

	for i in range(k):
		a = randrange(2, n - 1)
		if not check(a, s, d, n):
			return False
	return True

def is_prime(n):
    return miller_rabin(n, 30)

def get_prime(r0, k0):
    """Returns a prime within [r0*k0, 2**b] where b is the bit length of r0*k0"""
    prod = r0*k0
    b = int(log(prod, 2))
    while True:
        n = randrange(prod, 2*2**b)
        if is_prime(n):
            if debug:
                print("get_prime: return", n)
            return n

def get_parameters(p):
    """Needs improvement based on ECDSA paper"""
    # seed = randint(1, 2*p)
    a = randrange(0, p)
    b = randrange(0, p)
    return (a, b)

def schoof(a, b, p):
    """Ugly hack to call pyschoof and grab its output.
    Importing the file and calling pyschoof would be ideal,
    but hard to do with dependencies/path-related issues."""
    if debug:
        print("Schoof: points on", (a, b, p))
    result = subprocess.run(['python', './schoof.py', str(p), str(a), str(b)], stdout=subprocess.PIPE)
    output = str(result.stdout).split()[-1][:-3]
    if output == '':
        return 0
    return int(output)

def is_strong(r0, k0, p, N):
    """Given a curve defined over F_p of group order N,
    return 0 if it fits the security requirements of r0, k0
    return r such that |E(F_p)| = k*r"""
    if abs(N-(p+1)) > 2*sqrt(p):
        return 0
    r, k = 0, 0
    for i in range(1, k0+1):
        if N % i == 0 and is_prime(N/i) and N/i >= r0:
            r = N/i
            k = i
    if r == 0 or r == p:
        return 0
    p_r = 1
    for i in range(1, 20):
        p_r = p*p_r % r
        if p_r == 1:
            return 0
    return r

def random_approach(r0, k0):
    p = get_prime(r0, k0)
    while True:
        (a, b) = get_parameters(p)
        N = schoof(a, b, p)
        r = is_strong(r0, k0, p, N)
        if r != 0:
            return (a, b, p, r, N/r)

r0 = 2**10
k0 = 4

C = random_approach(r0, k0)
print(C)
