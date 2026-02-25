"""
    Module to compute the non-commuting charges of a DLA.
"""
from itertools import combinations
from paulie.common.pauli_string_collection import PauliStringCollection

def non_commuting_charges(generators: PauliStringCollection)->PauliStringCollection:
    """
    Finds the non-commuting charges of a DLA.

    Non-commuting charges describe non-Abelian symmetries — that is, elements of the
    stabilizer of the DLA that do not commute with each other.

    Args:
        generators (PauliStringCollection): Generating set of the DLA.
    Returns:
        PauliStringCollection: Set of non-commuting charges.
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
