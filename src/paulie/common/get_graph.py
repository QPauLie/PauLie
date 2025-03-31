def get_graph(generators, commutators=[]):
    """
    Get graph
    Args:
        generators: Array of PauliString
        commutators: The area of Pauli strings over which to build a graph. If not specified, then that's it
    Returns the vertices, edges, and labels of edges
    """
    vertices = []
    edge_labels = {}
    edges = []
    for a in generators:
        vertices.append(str(a))
        for b in generators:
            if a < b:
                if not a|b:
                    c = a^b
                    if len(commutators) == 0 or c in commutators:
                        edges.append((str(a), str(b)))
                        edge_labels[(str(a), str(b))] = str(c)
    return vertices, edges, edge_labels
