"""
Plot anti-commutation graph after transforming graph to a canonical type
"""
from paulie.helpers.drawing import plot_graph
from paulie.common.pauli_string_collection import PauliStringCollection


def plot_anti_commutation_graph(generators:PauliStringCollection) -> None:
    """
    Plot anti-commutation graph after tranform graph to canonic
    Args:
        generators: collection of Pauli strings
    Returns:
        None
    """
    vertices, edges, edge_labels =  generators.get_canonic_graph()
    plot_graph(vertices, edges, edge_labels)
