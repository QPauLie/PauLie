import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation

from paulie.helpers.recording import RecordGraph
from paulie.common.pauli import get_array_pauli_arrays
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


