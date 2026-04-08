"""
    Get the anticommutation graph for a set of generators.
"""
from paulie.common.pauli_string_bitarray import PauliString

def get_graph(generators: list[PauliString], commutators: list[PauliString] | None = None,
              flag_labels: bool = True) -> (tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]
                                            | tuple[list[str], list[tuple[str, str]]]):
    """
    Get the anticommutation graph.

    Args:
        generators (list[PauliString]): Generating set for anticommutator graph.
        commutators (list[PauliString], optional): Only add those edges which have labels belonging
            to this set. Defaults to None, in which case all edges are added.
        flag_labels (bool, optional): Boolean indicating whether to include labels in the graph.
            Defaults to `True`.
    Returns:
        Vertices, edges, and edge labels for the anticommutation graph.
    """
    if commutators is not None:
        commutators_set = set(commutators)
    else:
        commutators_set = None

    if hasattr(generators, 'generators'):
        generators = generators.generators

    vertices = [str(g) for g in generators]
    edge_labels = {}
    edges = []

    # Build overlap index
    qubit_to_indices: dict[int, list[int]] = {}
    for i, g in enumerate(generators):
        for qubit in g.get_support():
            if qubit not in qubit_to_indices:
                qubit_to_indices[qubit] = []
            qubit_to_indices[qubit].append(i)

    # Check only overlapping pairs
    for i, a in enumerate(generators):
        # Find all candidates that overlap with 'a'
        candidates = set()
        for qubit in a.get_support():
            if qubit in qubit_to_indices:
                for j in qubit_to_indices[qubit]:
                    if j > i:
                        candidates.add(j)

        for j in candidates:
            b = generators[j]
            c = a ^ b  # This returns None if they commute
            if c and (commutators_set is None or c in commutators_set):
                str_a, str_b = str(a), str(b)
                edges.append((str_a, str_b))
                if flag_labels:
                    edge_labels[(str_a, str_b)] = str(c)

    if flag_labels:
        return vertices, edges, edge_labels
    else:
        return vertices, edges
