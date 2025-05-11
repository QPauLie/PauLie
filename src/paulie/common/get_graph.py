"""
Get anticommutator graph
"""
from itertools import combinations

def get_graph(generators, commutators=None, flag_labels = True):
    """
    Get anticommutator graph
    Args:
        generators: Array of PauliString
        commutators: The area of Pauli strings over which to build a graph. 
        If not specified, then that's it
    Returns the vertices, edges, and labels of edges
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
