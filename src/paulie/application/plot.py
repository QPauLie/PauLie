"""
    Module to plot the anti-commutation graph after transforming it to a canonical type.
"""
from paulie.helpers.drawing import plot_graph
from paulie.common.pauli_string_collection import PauliStringCollection


def plot_anti_commutation_graph(generators:PauliStringCollection) -> None:
    """
    Plot a canonical anti-commutation graph of a set of generators.

    Args:
        generators (PauliStringCollection): Generating set of the Pauli string DLA.
    Returns:
        None
    """
    vertices, edges, edge_labels =  generators.get_canonic_graph()
    plot_graph(vertices, edges, edge_labels)
