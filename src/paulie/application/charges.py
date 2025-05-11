"""
    Non-commuting charges describe non-Abelian symmetries — that is, elements of the
    stabilizer of the DLA that do not commute with each other.
"""
from itertools import combinations
from paulie.common.pauli_string_collection import PauliStringCollection

def non_commuting_charges(generators: PauliStringCollection)->PauliStringCollection:
    """
    Non-commuting charges describe non-Abelian symmetries — that is, elements of the
    stabilizer of the DLA that do not commute with each other.

    inputs: generators as strings
    outputs: list of charges as strings
    """
    non_q = PauliStringCollection()
    comm = generators.get_commutants()
    for c,q in combinations(comm,2):
        if c | q is False:
            if c not in non_q:
                non_q.append(c)
            if q not in non_q:
                non_q.append(q)
    return non_q
