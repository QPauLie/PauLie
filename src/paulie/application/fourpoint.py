"""
    Module to compute the four-point correlator of Pauli strings.
"""
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
from paulie.application.otoc import average_otoc


def fourpoint(generators: PauliStringCollection,
            p: PauliString,
            q: PauliString,
            r: PauliString,
            s: PauliString) -> float:
    r"""
    Computes the four-point correlator of Pauli strings :math:`P`, :math:`Q`, 
    :math:`R`, and :math:`S`.
    
    For :math:`PR,QS \propto L` where :math:`L` is a linear symmetry, this
    reduces to the average OTOC of :math:`P` and :math:`Q`.

    (arXiV: 2502.16404)

    Args:
        generators (PauliStringCollection): Generating set of the Pauli string DLA.
        p (PauliString): Pauli string :math:`P`.
        q (PauliString): Pauli string :math:`Q`.
        r (PauliString): Pauli string :math:`R`.
        s (PauliString): Pauli string :math:`S`.
    Returns:
        float: The four-point correlator of the Pauli strings.
    """
    commutant = generators.get_commutants()
    rp = r@p
    qs = q@s
    if rp == qs and qs in commutant:
        return average_otoc(generators, p, q)
    return 0
