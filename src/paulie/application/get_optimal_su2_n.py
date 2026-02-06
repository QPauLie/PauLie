"""Get generator set with optimal generation rate for su(2**n)"""
import math
from paulie.common.pauli_string_collection import PauliStringCollection


def get_optimal_connections_su_2_n(generators:PauliStringCollection) -> tuple[int, int]:
    r""" Get optimal number of connections in anticommutation graph for
    :math:`\mathfrak{su}(2^{n})`."""
    fraction = 0.706
    ng = len(generators)
    total_pair = ng*(ng-1)/2
    if total_pair < 1:
        return -1
    n_pair = fraction * total_pair
    n_pair = math.floor(n_pair)
    return n_pair

def get_optimal_su_2_n_generators(generators:PauliStringCollection) -> PauliStringCollection|None:
    r""" Get optimal generator set for :math:`\mathfrak{su}(2^{n})`."""
    g = generators.copy().get_independents()
    n_pair = get_optimal_connections_su_2_n(g)
    return g.find_generators_with_connection(n_pair) if n_pair > -1 else None
