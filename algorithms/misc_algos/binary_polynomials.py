#!/bin/python

from numpy import polymul

# For polynomials in F2[X]/(X^3+X+1)

X3 = 2**3
X2 = 2**2
X = 2

P = X3 + X + 1

def repr_poly(p):
    return "{0:b}".format(p)

def str_poly(p):
    if p == 0:
        return '0'

    bstr =  repr_poly(p)
    s = ''
    for i, c in enumerate(bstr[::-1]):
        if c == "0":
            continue
        if i == 0:
            s += "1"
        if i == 1:
            s += "X"
        if i >= 2:
            s += "X^"+str(i)
        if i != len(bstr) - 1:
            s += " + "
    return s

def deg_poly(p):
    return len(repr_poly(p)) - 1

def add_poly(p1, p2):
    return p1 ^ p2

def sub_poly(p1, p2):
    return p1 ^ p2

def poly_to_array(p):
    bstr = repr_poly(p)
    return [int(c) for c in bstr]

def array_to_poly(arr):
    p = 0
    for i, c in enumerate(arr[::-1]):
        p += c*2**i
    return p

def mul_poly(p1, p2):
    p3 = polymul(poly_to_array(p1), poly_to_array(p2))
    p3 = [c % 2 for c in p3]
    p3 = array_to_poly(p3)
    while deg_poly(p3) >= 3:
        shift = deg_poly(p3) - 3
        p3 = p3 ^ (P << shift)
    return p3

print(str_poly(mul_poly(1+X+X2, 1+X+X2)))
