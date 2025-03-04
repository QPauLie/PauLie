from paulie.helpers.classify_generators import classify_generators
from paulie.graphs.graph_view import get_graph_view
from paulie.helpers.drawing import plot_graph


# Plot anti-commutation graph after tranform graph to canonic
# generators - list of generators
# size - Generator extensions to size size
def plot_anti_commutation_graph(generators, size=0):
    classification = classify_generators(generators, size=size)
    vertices, edges, edge_labels =  get_graph_view(classification.get_vertices())
    plot_graph(vertices, edges, edge_labels)