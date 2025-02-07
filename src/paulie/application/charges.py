from paulie.common.pauli import commutant, is_commutative, get_array_pauli_arrays,get_array_pauli_strings
from itertools import combinations


def non_commuting_charges(generators):
    """
    inputs: generators as strings
    outputs: list of charges as strings
    """
    non_q = []
    comm = get_array_pauli_arrays(commutant(generators))
    for c,q in combinations(comm,2):
        if is_commutative(c, q) is False:
            if c not in non_q:
                non_q.append(c)
            if q not in non_q:
                non_q.append(q)
    return get_array_pauli_strings( non_q)
