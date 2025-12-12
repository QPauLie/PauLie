"""
Get anticommutator graph
"""
from itertools import combinations
from paulie.common.pauli_string_bitarray import PauliString

def get_graph(generators:list[PauliString], commutators:list[PauliString]|None=None
, flag_labels:bool = True
) -> (tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]
     |tuple[list[str], list[tuple[str, str]]]):
    """
    Get anticommutator graph.

    Args:
        generators: Array of PauliString.
        commutators: Area of Pauli strings over which to build a graph.
        flag_labels: Boolean indicating whether to include labels in the graph.
    Returns:
        Vertices, edges, and labels of edges.
    """
    if not commutators:
        commutators = []
    vertices = [str(g) for g in generators]
    edge_labels = {}
    edges = []
    for a, b in combinations(generators, 2):
        c = a^b
        if c and (len(commutators) == 0 or c in commutators):
            edges.append((str(a), str(b)))
            if flag_labels:
                edge_labels[(str(a), str(b))] = str(c)
    if flag_labels:
        return vertices, edges, edge_labels
    else:
        return vertices, edges
