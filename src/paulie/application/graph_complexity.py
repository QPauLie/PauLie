"""
    Computes the average graph complexity of a specified PauliString.
"""
import networkx as nx
from paulie.common.pauli_string_collection import PauliString, PauliStringCollection

def average_graph_complexity(generators: PauliStringCollection, p: PauliString):
    """
    Return the average graph complexity 
    """
    # Get commutator graph
    vertices, edges = generators.get_commutator_graph()
    # Construct graph in NetworkX
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    graph.add_edges_from(edges)
    # Get connected component containing p
    subgraph = graph.subgraph(nx.node_connected_component(graph, str(p)))
    # Return a dict of shortest path lengths to all nodes in connected component from p
    spl = nx.shortest_path_length(subgraph, source=str(p))
    # Sum the path lengths and divide by connected component size to obtain average graph complexity
    return sum(spl.values()) / subgraph.number_of_nodes()
