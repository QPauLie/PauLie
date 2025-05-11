"""
Plot anti-commutation graph after tranform graph to canonic
"""
from paulie.helpers.drawing import plot_graph
from paulie.common.pauli_string_collection import PauliStringCollection


def plot_anti_commutation_graph(generators:PauliStringCollection):
    """
    Plot anti-commutation graph after tranform graph to canonic
    """
    vertices, edges, edge_labels =  generators.get_canonic_graph()
    plot_graph(vertices, edges, edge_labels)
