#!/bin/env python

from paulie.common.pauli_string_factory import get_pauli_string as ps
from paulie.application.otoc import average_otoc
from paulie.application.fourpoint import fourpoint

# These are the generators of n=3 "matchgate" dynamics
generators = ps(["ZII", "IZI", "IIZ", "XXI", "IXX"])
p = ps("XYZ")
q = ps("YZZ")
r = ps("YXZ")
s = ps("YYZ")

print(f"G:{generators}")
print(f"P:{p}")
print(f"Q:{q}")
print(f"R:{r}")
print(f"S:{s}")
print(f"Average OTOC of PQ: {average_otoc(generators,p,q)}")
print(f"Fourpoint for R=P and S=Q: {fourpoint(generators,p,q,p,q)}")
print(f"Fourpoint for PQRS: {fourpoint(generators,p,q,r,s)}")
