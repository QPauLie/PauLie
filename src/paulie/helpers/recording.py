from paulie.common.get_graph import get_graph

class FrameGraph:
      def __init__(self, vertices, edges, edge_labels = None):
          self.vertices = vertices
          self.edges = edges
          self.edge_labels = edge_labels

      def get_graph(self):
          return self.vertices.copy(), self.edges.copy(), self.edge_labels

class FrameRecord:
      def __init__(self, graph: FrameGraph=None, lighting=None, appending=None, contracting=None, lits=None, p=None, q=None, removing_vertices=None, replacing_vertices=None, dependent=False, title=None, init=False):
          self.graph = graph
          self.lighting = str(lighting) if lighting else None
          self.lits = [str(v) for v in lits] if lits else None
          self.p = str(p) if p else None
          self.q = str(q) if q else None
          self.removing_vertices = [str(v) for v in removing_vertices] if removing_vertices else None
          self.replacing_vertices = [str(v) for v in replacing_vertices] if replacing_vertices else None
          self.dependent = str(dependent) if dependent else None
          self.title = title
          self.contracting = str(contracting) if contracting else None
          self.appending = str(appending) if appending else None
          self.init = init


      def get_graph(self):
          if not self.graph:
              return None
          return self.graph.get_graph()

      def get_lighting(self):
          return self.lighting

      def get_title(self):
          return self.title

      def is_appending(self):
          return not self.appending

      def get_is_appending(self, vertix):
          if not self.appending:
              return False
          return self.appending == vertix

      def get_is_contracting(self, vertix):
          if not self.contracting:
              return False
          return self.contracting == vertix

      def get_is_p(self, vertix):
          if not self.p:
              return False
          return self.p == vertix

      def get_is_q(self, vertix):
          if not self.q:
              return False
          return self.q == vertix

      def get_is_dependent(self, vertix):
          if not self.dependent:
              return False
          return self.dependent == vertix

      def get_is_lits(self, vertix):
          if not self.lits:
              return False
          return vertix in self.lits

      def is_removing(self):
          return not self.removing_vertices

      def get_is_removing(self, vertix):
          if not self.removing_vertices:
              return False
          return vertix in self.removing_vertices

      def get_is_replacing(self, vertix):
          if not self.replacing_vertices:
              return False
          return vertix in self.replacing_vertices

      def get_init(self):
          return self.init




class RecordGraph: 
     def __init__(self):
         self.frames: FrameRecord = []
         self.positions = {}
         self.x_position_lighting = 0


     def append_frame(self, frame: FrameRecord):
         self.frames.append(frame)

     def append(self, graph: FrameGraph=None, lighting=None, appending=None, contracting=None, lits=None, p=None, q=None, removing_vertices=None, replacing_vertices=None, dependent=False, title=None, init=False):
         self.append_frame(FrameRecord(graph, lighting=lighting, appending=appending, contracting=contracting, lits=lits, p=p, q=q, removing_vertices=removing_vertices, replacing_vertices=replacing_vertices, dependent=dependent, title=title, init=init))


     def get_frame(self, index):
         if index > len(self.frames) - 1:
            raise ValueError("Out of index")
         return self.frames[index]

     def clear(self):
         self.frames = []

     def get_size(self):
         return len(self.frames)

     def get_graph(self, index):
         while index > -1:
             frame = self.get_frame(index)
             graph = frame.get_graph()
             if not graph:
                 index -= 1
                 continue
             return graph
         return None

     def get_is_prev(self, index):
         frame = self.get_frame(index)
         return frame.get_graph() is None

     def set_positions(self, positions):
         self.positions = positions

     def get_positions(self):
         return self.positions

     def set_x_position_lighting(self, x_position_lighting):
         self.x_position_lighting = x_position_lighting

     def get_x_position_lighting(self):
         return self.x_position_lighting


def recording_graph(record:RecordGraph, collection=None, lighting=None, appending=None, contracting=None, lits=None, p=None, q = None, removing_vertices=None, replacing_vertices=None, dependent=False, title=None, init=False):

    graph = None
    #vertices=[]
    #edges = []
    #edge_labels = []
    if collection is not None:
        vertices, edges, edge_labels = get_graph(collection)
        graph = FrameGraph(vertices, edges, edge_labels)
    record.append(graph, lighting=lighting, appending=appending, contracting=contracting, lits=lits, p=p, q=q, removing_vertices=removing_vertices, replacing_vertices=replacing_vertices, dependent=dependent, title=title, init=init)


