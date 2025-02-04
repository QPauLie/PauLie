from paulie.helpers.transform import app_transform_to_canonics
from paulie.classifier.transform import merge_canonics
from paulie.graphs.graph_view import get_graph_view
from paulie.helpers.drawing import plot_graph


# Plot anti-commutation graph after tranform graph to canonic
# generators - list of generators
# size - Generator extensions to size size
def plot_anti_commutation_graph(generators, size=0):
    canonics = app_transform_to_canonics(generators, size=size)
    nodes = merge_canonics(canonics)
    vertices, edges, edge_labels =  get_graph_view(nodes)
    plot_graph(vertices, edges, edge_labels)