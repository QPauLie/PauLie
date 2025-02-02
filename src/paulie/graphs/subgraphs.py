import networkx as nx

from paulie.common.pauli import get_pauli_array
from paulie.graphs.graph_view import get_graph_view


def get_subgraphs(nodes, commutators = []):
    subgraphs = []
    vertices, edges, _ = get_graph_view(nodes, commutators)

    g = nx.Graph()
    g.add_nodes_from(vertices)
    g.add_edges_from(edges)
    string_subgraphs = sorted(nx.connected_components(g), key=len, reverse=True)

    for sb in string_subgraphs:
        subgraph = [] 
        for item in sb:
            subgraph.append(get_pauli_array(item))
        subgraphs.append(subgraph)
    return subgraphs
