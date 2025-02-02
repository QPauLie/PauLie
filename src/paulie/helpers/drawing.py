import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation

from paulie.helpers.recording import RecordGraph, recording_graph
from paulie.common.pauli import get_array_pauli_arrays
from paulie.classifier.transform import merge_canonics, transform_to_canonics
from paulie.common.ext_k_local import get_k_local_generators
from paulie.graphs.graph_view import get_graph_view


def plot_graph(vertices, edges, edge_labels = None):
   graph = nx.Graph()
   graph.add_nodes_from(vertices)
   graph.add_edges_from(edges)
   pos = nx.spring_layout(graph)
   if edge_labels is not None:
       nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=edge_labels,font_color="red")
   nx.draw_networkx(graph, pos=pos)
   plt.show()


def plot_graph_by_nodes(nodes, commutators=[]):
    vertices, edges, edge_labels = get_graph_view(nodes, commutators)
    return plot_graph(vertices, edges, edge_labels)


def plot_graph_by_string_nodes(nodes, commutators=[]):
    vertices, edges, edge_labels = get_graph_view(get_array_pauli_arrays(nodes), get_array_pauli_arrays(commutators))
    return plot_graph(vertices, edges, edge_labels)


def animation_graph(record: RecordGraph, interval=1000, repeat=False, storage=None):
   graph = nx.Graph()
   fig, ax = plt.subplots(figsize=(6,4))
   def clear():
       ax.clear()
       graph.remove_nodes_from(list(n for n in graph.nodes)) 

   def update(num):
       clear()
       vertices, edges, edge_labels = record.get_frame(num).get_graph()
       graph.add_nodes_from(vertices)
       graph.add_edges_from(edges)
       pos = nx.spring_layout(graph)
       if edge_labels is not None:
           nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=edge_labels,font_color='red')
       return nx.draw_networkx(graph, pos=pos)


   ani = matplotlib.animation.FuncAnimation(fig, update, frames=record.get_size(), interval=interval, repeat=repeat)
   if storage is not None:
       ani.save(filename=storage["filename"], writer=storage["writer"])
   plt.show()


#Partial function for constructing a canonical graph
# generators - list of generators
# size - Generator extensions to size size
# debug - debbuging
# record - recording 
def _transform_to_canonics(generators, size = 0, debug=False, record=None, initGraph=False):
    if size == 0:
        bitGenerators = get_array_pauli_arrays(generators)
    else:
        bitGenerators = get_k_local_generators(size, generators)
    if record is not None and initGraph:
        recording_graph(record, bitGenerators)
    return transform_to_canonics(bitGenerators, debug, record)

# Get algebra
# generators - list of generators
# size - Generator extensions to size size
def get_algebra(generators, size=0):
    canonics = _transform_to_canonics(generators, size=size)
    algebras = []
    for canonic in canonics:
        algebras.append(canonic["shape"].get_algebra())
    return " + ".join(algebras)


# Animation building transformation anti-commutation graph
# generators - list of generators
# size - Generator extensions to size size
def animation_anti_commutation_graph(generators, size=0, storage=None, interval=1000, initGraph=False):
    record = RecordGraph()
    animation_graph(record, storage=storage, interval=interval)


# Plot anti-commutation graph after tranform graph to canonic
# generators - list of generators
# size - Generator extensions to size size
def plot_anti_commutation_graph(generators, size=0):
    canonics = _transform_to_canonics(generators, size=size)
    nodes = merge_canonics(canonics)
    vertices, edges, edge_labels =  get_graph_view(nodes)
    plot_graph(vertices, edges, edge_labels)