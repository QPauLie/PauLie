from paulie.common.ext_k_local import (
    get_k_local_algebra_generators,
    get_k_local_nested_nodes_in_algebra
)
from paulie.helpers.drawing import plot_graph_by_nodes


if __name__ == "__main__":
    generators = get_k_local_algebra_generators(2, "a13")
    nestedNodes = get_k_local_nested_nodes_in_algebra(2, "a13")
    plot_graph_by_nodes(nestedNodes, generators)
