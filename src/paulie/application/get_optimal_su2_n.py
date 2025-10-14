"""Get optimal generator set by fraction for su(2**n)"""
import math
from paulie.common.pauli_string_collection import PauliStringCollection


def get_optimal_connections_su_2_n(generators:PauliStringCollection, fraction: float
    )->tuple[int, int]:
    """Get optimal connections in graph for su(2**n)"""
    ng = len(generators)
    total_pair = ng*(ng-1)/2
    if total_pair < 1:
        return (-1, -1)
    n_pair = fraction * total_pair
    n1_pair = math.floor(n_pair)
    n2_pair = math.floor(n_pair + 1)
    return (n1_pair, n2_pair)

def get_optimal_su_2_n_generators_by_fraction_left(generators:PauliStringCollection, fraction: float
    ) -> PauliStringCollection|None:
    """ Get optimal generator set by fraction for su(2**n) round left"""
    g = generators.copy().get_independents()
    n_pair = get_optimal_connections_su_2_n(g, fraction)[0]
    return g.find_generators_with_connection(n_pair) if n_pair > -1 else None

def get_optimal_su_2_n_generators_by_fraction_right(generators:PauliStringCollection,
    fraction: float) -> PauliStringCollection|None:
    """ Get optimal generator set by fraction for su(2**n) round right"""
    g = generators.copy().get_independents()
    n_pair = get_optimal_connections_su_2_n(g, fraction)[1]
    return g.find_generators_with_connection(n_pair) if n_pair > -1 else None
