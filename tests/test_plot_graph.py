"""Test plot_graph_by_nodes"""
import matplotlib.pyplot as plt
from paulie import (
    plot_graph_by_nodes,
    get_pauli_string as p,
    plot_anti_commutation_graph
)

def test_plot_graph() -> None:
    """Test plotting"""
    generators = p(["XIIII", "ZIIII", "IYXII", "IXIXI", "IIIZI", "IZIXZ", "IIIIX"])
    with plt.ion():
        plot_graph_by_nodes(generators)
        assert True
    with plt.ion():
        plot_anti_commutation_graph(generators)
        assert True
