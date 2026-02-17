"""
Plot anti-commutation graph after transforming graph to a canonical type.
"""
from paulie.helpers.drawing import plot_graph
from paulie.common.pauli_string_collection import PauliStringCollection


def plot_anti_commutation_graph(generators:PauliStringCollection) -> None:
    """
    Plot a canonical anti-commutation graph of the passed generators.

    Args:
        generators: Collection of Pauli strings.
    Returns:
        None
    """
    vertices, edges, edge_labels =  generators.get_canonic_graph()
    plot_graph(vertices, edges, edge_labels)
