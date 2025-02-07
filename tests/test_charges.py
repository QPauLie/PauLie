from paulie.common.algebras import get_lie_algebra
from paulie.application.charges import non_commuting_charges


# a_8, a_9 posses non-commuting charges for all n
generators = get_lie_algebra()["a8"]
assert len(non_commuting_charges(generators)) != 0
generators = get_lie_algebra()["a9"]
assert len(non_commuting_charges(generators)) != 0
