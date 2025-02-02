from paulie.common.pauli import (
    get_array_pauli_arrays, 
    get_pauli_string,
    is_commutate,
    multi_pauli_arrays,
)


def get_graph_view(nodes, commutators=[]):
    vertices = []
    edge_labels = {}
    edges = []
    for nodeA in nodes:
        vertices.append(get_pauli_string(nodeA))
        for nodeB in nodes:
            if nodeA < nodeB:
                if is_commutate(nodeA, nodeB) is False:
                    nodeC = multi_pauli_arrays(nodeA, nodeB)
                    if len(commutators) == 0 or nodeC in commutators:
                        edges.append((get_pauli_string(nodeA), get_pauli_string(nodeB)))
                        edge_labels[(get_pauli_string(nodeA), get_pauli_string(nodeB))] = get_pauli_string(nodeC)
    return vertices, edges, edge_labels


def get_graph_view_by_string(nodes, commutators=[]):
    return get_graph_view(get_array_pauli_arrays(nodes), get_array_pauli_arrays(commutators))
