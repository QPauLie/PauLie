"""Test plot_graph_by_nodes"""
import matplotlib.pyplot as plt
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.helpers.drawing import plot_graph_by_nodes


def test_plot_graph() -> None:
    """Test subgraph"""
    generators = p(["XIIII", "ZIIII", "IYXII", "IXIXI", "IIIZI", "IZIXZ", "IIIIX"])
    with plt.ion():
        plot_graph_by_nodes(generators)
        assert True
