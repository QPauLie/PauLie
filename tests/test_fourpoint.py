"""
    Tests for four point.
"""

from paulie.common.pauli_string_factory import get_pauli_string as ps
from paulie.application.otoc import average_otoc
from paulie.application.fourpoint import fourpoint

# These are the generators of n=3 "matchgate" dynamics
generators = ps(["ZII", "IZI", "IIZ", "XXI", "IXX"])
p = ps("XYZ")
q = ps("YZZ")
r = ps("YXZ")
s = ps("YYZ")

assert fourpoint(generators,p,q,p,q) == average_otoc(generators,p,q)
assert fourpoint(generators,p,q,r,s) == 0
