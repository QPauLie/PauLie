"""Test plot_graph_by_nodes"""
import matplotlib.pyplot as plt
from paulie import plot_graph_by_nodes, get_pauli_string as p


def test_plot_graph() -> None:
    """Test subgraph"""
    generators = p(["XIIII", "ZIIII", "IYXII", "IXIXI", "IIIZI", "IZIXZ", "IIIIX"])
    with plt.ion():
        plot_graph_by_nodes(generators)
        assert True
