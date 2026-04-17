"""
    Module with graph drawing utilities.
"""
import math
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.get_graph import get_graph

def plot_graph(vertices:list[str],
               edges:list[tuple[str,str]],
               edge_labels:dict[tuple[str,str],str] = None) -> None:
    """
    Plot an arbitrary graph.

    Args:
        vertices (list[str]): List of vertices.
        edges (list[tuple[str, str]]): List of edges.
        edge_labels (dict[tuple[str,str],str], optional): List of edge labels. Defaults to None, in
            which case no edge labels are drawn.
    Returns:
        None
    """
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    graph.add_edges_from(edges)
    pos = nx.spring_layout(graph)
    if edge_labels is not None:
        nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=edge_labels,font_color="red")
    nx.draw_networkx(graph, pos=pos)
    #plt.show()

def plot_graph_by_nodes(nodes:PauliStringCollection,
                        commutators:PauliStringCollection|list[PauliString]=None
) -> None:
    """
    Plot anticommutator graph with edges labeled by commutator of endpoints.

    Args:
        nodes (PauliStringCollection): Generating set for anticommutator graph.
        commutators (PauliStringCollection|list[PauliString], optional): Only show those edges which
            have labels belonging to this set. Defaults to None, in which case all edges are shown.
    Returns:
        None
    """
    if not commutators:
        commutators = []
    vertices, edges, edge_labels = get_graph(nodes, commutators)
    return plot_graph(vertices, edges, edge_labels)
