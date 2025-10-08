"""
    Compute the average out-of-time-order correlator between two Pauli strings.
"""
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
from paulie.application.otoc import average_otoc


def fourpoint(generators: PauliStringCollection,
            p: PauliString,
            q: PauliString,
            r: PauliString,
            s: PauliString) -> float:
    """
    Computes the four point according to Proposition 2. 
    For PQ, RS prop L this reduces to the average OTOC of P and Q.

    (arXiV: 2502.16404)

    Args:
        generators: Generating set of the Pauli string DLA.
        p: Pauli string P
        q: Pauli string Q
        r: Pauli string R
        s: Pauli string S
    """
    l = generators.get_commutants()
    rp = r@p
    qs = q@s
    if rp == qs and qs in l:
        return average_otoc(generators, p, q)
    return 0
