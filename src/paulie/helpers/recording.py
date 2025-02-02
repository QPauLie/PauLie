

from paulie.common.pauli import get_array_pauli_arrays
from paulie.graphs.graph_view import get_graph_view


class FrameGraph:
      def __init__(self, vertices, edges, edge_labels = None):
          self.vertices = vertices
          self.edges = edges
          self.edge_labels = edge_labels

      def get_graph(self):
          return self.vertices, self.edges, self.edge_labels


class RecordGraph: 
     def __init__(self):
         self.frames: FrameGraph = []

     def append_frame(self, frame: FrameGraph):
         self.frames.append(frame)

     def append(self, vertices, edges, edge_labels = None):
         self.append_frame(FrameGraph(vertices, edges, edge_labels))

     def get_frame(self, index):
         if index > len(self.frames) - 1:
            return FrameGraph(["None"], None)
         return self.frames[index]

     def clear(self):
         self.frames = []

     def get_size(self):
         return len(self.frames)


def recording_graph(record: RecordGraph, nodes):
    if record is None:
        return
    vertices, edges, edge_labels = get_graph_view(nodes)
    record.append(vertices, edges, edge_labels)


def recording_graph_string(record: RecordGraph, nodes):
    if record is None:
        return
    vertices, edges, edge_labels = get_graph_view(get_array_pauli_arrays(nodes))
    record.append(vertices, edges, edge_labels)
