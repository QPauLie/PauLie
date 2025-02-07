from paulie.common.algebras import G_LIE
from paulie.application.charges import non_commuting_charges


# a_8, a_9 posses non-commuting charges for all n
generators = G_LIE["a8"]
assert len(non_commuting_charges(generators)) != 0
generators = G_LIE["a9"]
assert len(non_commuting_charges(generators)) != 0
