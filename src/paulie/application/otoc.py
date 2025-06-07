"""
    Compute the average out-of-time-order correlator between two Pauli strings.
"""
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
import networkx as nx

def average_otoc(generators: PauliStringCollection, v: PauliString, w: PauliString):
    # Generate commutator graph
    vertices, edges = generators.get_commutator_graph()
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    graph.add_edges_from(edges)
    # Get connected component of V
    v_connected_component = nx.node_connected_component(graph, str(v))
    # Count the number of elements t in the connected component of V
    # that anticommute with W
    anticommute_count = sum([not w | p(t) for t in v_connected_component])
    return 1 - 2 * anticommute_count / len(v_connected_component)
