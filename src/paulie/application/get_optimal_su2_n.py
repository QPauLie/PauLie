"""
    Module to get the generator set with optimal generation rate for su(2**n)
"""
import math
from paulie.common.pauli_string_collection import PauliStringCollection


def get_optimal_edges_su_2_n(ng:int) -> int:
    r"""
    Get optimal number of edges in anticommutation graph for :math:`\mathfrak{su}(2^{n})`.

    Args:
        ng (int): Number of generators.
    Returns:
        int: Optimal number of edges in anticommutation graph for :math:`\mathfrak{su}(2^{n})`.
    """
    fraction = 0.706
    total_pair = ng*(ng-1)/2
    if total_pair < 1:
        return -1
    n_pair = fraction * total_pair
    n_pair = math.floor(n_pair)
    return n_pair

def get_optimal_su_2_n_generators(generators:PauliStringCollection) -> PauliStringCollection|None:
    """ Get optimal generator set for :math:`\\mathfrak{su}(2^{n})`."""
    g = generators.copy().get_independents()
    n_pair = get_optimal_edges_su_2_n(len(g))
    return g.find_generators_with_connection(n_pair) if n_pair > -1 else None
