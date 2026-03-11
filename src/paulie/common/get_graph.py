"""
    Get the anticommutation graph for a set of generators.
"""
from itertools import combinations
from paulie.common.pauli_string_bitarray import PauliString

def get_graph(generators:list[PauliString], commutators:list[PauliString]=None
, flag_labels:bool=True
) -> (tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]
     |tuple[list[str], list[tuple[str, str]]]):
    """
    Get the anticommutation graph.

    Args:
        nodes (list[PauliString]): Generating set for anticommutator graph.
        commutators (list[PauliString], optional): Only add those edges which have labels belonging
            to this set. Defaults to None, in which case all edges are added.
        flag_labels (bool, optional): Boolean indicating whether to include labels in the graph.
            Defaults to `True`.
    Returns:
        Vertices, edges, and edge labels for the anticommutation graph.
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
