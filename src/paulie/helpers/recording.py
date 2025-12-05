"""
Recording of building graph
"""
import numpy as np
from paulie.common.get_graph import get_graph
from paulie.common.pauli_string_bitarray import PauliString

class FrameGraph:
    """
    Frame
    """
    def __init__(self, vertices:list[str], edges:list[tuple[str,str]],
                 edge_labels:dict[tuple[str,str],str] = None
    ) -> None:
        """
        Constructor.

        Args:
            vertices: List of vertices.
            edges: List of edges.
            edge_labels: List of edge's labels.
        Returns:
            None
        """
        self.vertices = vertices
        self.edges = edges
        self.edge_labels = edge_labels

    def get_graph(self
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get graph pf the frame.

        Returns:
            Vertices, edges, and labels of edges.
        """
        return self.vertices.copy(), self.edges.copy(), self.edge_labels

class FrameRecord:
    """
    Frame of record.
    """
    def __init__(self, graph:FrameGraph=None, lighting:PauliString=None,
                 appending:PauliString=None,  contracting=PauliString,
                 lits:list[PauliString]=None, p:PauliString=None,
                 q:PauliString=None, removing_vertices:PauliString=None,
                 replacing_vertices:list[PauliString]=None, dependent:PauliString=None,
                 title:str=None, init:bool=False) -> None:
        """
        Constructor.

        Args:
            graph: Frame of graph.
            lighting: Canonical graph join candidate.
            appending: Pauli string for appending.
            contracting: Pauli string for contracting.
            lits: List of lited vertices.
            p: Lited Pauli string in legs long 1.
            q: Unlited Pauli string in legs long 1.
            removing_vertices: List of removing vertices.
            replacing_vertices: List of replacing vertices.
            dependent: Dependent vertex.
            title: Title of frame.
            init: Init state of building.
        Returns:
            None
        """
        self.graph = graph
        self.lighting = str(lighting) if lighting else None
        self.lits = [str(v) for v in lits] if lits else []
        self.p = str(p) if p else None
        self.q = str(q) if q else None
        self.removing_vertices = [str(v) for v in removing_vertices] if removing_vertices else []
        self.replacing_vertices = ([str(v) for v in replacing_vertices]
                                   if replacing_vertices else [])
        self.dependent = str(dependent) if dependent else None
        self.title = title
        self.contracting = str(contracting) if contracting else None
        self.appending = str(appending) if appending else None
        self.init = init

    def get_graph(self
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get graph of record's frame.

        Returns:
            Vertices, edges, and labels of edges.
        """
        if not self.graph:
            return None
        return self.graph.get_graph()
    def get_lighting(self) -> str:
        """
        Get lighting.

        Returns:
            Lighting.
        """
        return self.lighting

    def get_title(self) -> str:
        """
        Get title.

        Returns:
            Title.
        """
        return self.title

    def is_appending(self) -> bool:
        """
        Check appending.

        Returns:
            True if not appending.
        """
        return not self.appending

    def get_is_appending(self, vertex: str) -> bool:
        """
        Check vertex is appending.

        Args:
            vertex: Vertex for checking.
        Returns:
            True if vertex is appending.
        """
        if not self.appending:
            return False
        return self.appending == vertex

    def get_is_contracting(self, vertex:str) -> bool:
        """
        Check vertex is contracting.

        Args:
            vertex: Vertex for checking.
        Returns:
            True if vertex is contracting.

        """
        if not self.contracting:
            return False
        return self.contracting == vertex

    def get_is_p(self, vertex:str) -> bool:
        """
        Check vertex is p.

        Args:
            vertex: vVrtex for checking.
        Returns:
            True if vertex is lited in leg long 1.
        """
        if not self.p:
            return False
        return self.p == vertex

    def get_is_q(self, vertex:str) -> bool:
        """
        Check vertex is q.

        Args:
            vertex: Vertex for checking.
        Returns:
            True if vertex is unlited in leg long 1.
        """
        if not self.q:
            return False
        return self.q == vertex

    def get_is_dependent(self, vertex:str) -> bool:
        """
        Check vertex is dependent.

        Args:
            vertex: Vertex for checking.
        Returns:
            True if vertex is dependent.
        """
        if not self.dependent:
            return False
        return self.dependent == vertex

    def get_is_lits(self, vertex:str) -> bool:
        """
        Check vertex is lited.

        Args:
            vertex: Vertex for checking.
        Returns:
            True if vertex is lited.
        """
        return vertex in self.lits

    def is_removing(self) -> bool:
        """
        Check is removing.

        Returns:
            True if vertex is not removing.
        """
        return not self.removing_vertices

    def get_is_removing(self, vertex:str) -> bool:
        """
        Check vertex is removing.

        Args:
            vertex: Vertex for checking.
        Returns:
            True if vertex is removing.
        """
        return vertex in self.removing_vertices

    def get_is_replacing(self, vertex:str) -> bool:
        """
        Check vertex is replacing.

        Args:
            vertex: Vertex for checking.
        Returns:
            True if vertex is replacing.
        """
        return vertex in self.replacing_vertices

    def get_init(self) -> bool:
        """
        Get init frame.

        Returns:
            True if init state.
        """
        return self.init

class RecordGraph:
    """
    Record of graph.
    """
    def __init__(self) -> None:
        """
        Constructor.

        Returns:
            None
        """
        self.frames: FrameRecord = []
        self.positions = {}
        self.x_position_lighting = 0


    def append_frame(self, frame: FrameRecord) -> None:
        """
        Append frame.

        Args:
            frame: Frame of record.
        Returns:
            None
        """
        self.frames.append(frame)

    def append(self, graph: FrameGraph=None,
               lighting:PauliString=None, appending:PauliString=None,
               contracting:PauliString=None, lits:list[PauliString]=None,
               p:PauliString=None, q:PauliString=None,
               removing_vertices:list[PauliString]=None, replacing_vertices:list[PauliString]=None,
               dependent:PauliString=None, title:str=None,
               init:bool=False) -> None:
        """
        Make frame and append to record.

        Args:
            graph: Frame of graph.
            lighting: Canonical graph join candidate.
            appending: Pauli string for appending.
            contracting: Pauli string for contracting.
            lits: List of lited vertices.
            p: Lited Pauli string in legs long 1.
            q: Unlited Pauli string in legs long 1.
            removing_vertices: List of removing vertices.
            replacing_vertices: List of replacing vertices.
            dependent: Dependent vertex.
            title: Title of frame.
            init: Init state of building.
        Returns:
            None
        """
        self.append_frame(FrameRecord(graph, lighting=lighting, appending=appending,
                                      contracting=contracting, lits=lits, p=p, q=q,
                                      removing_vertices=removing_vertices,
                                      replacing_vertices=replacing_vertices,
                                      dependent=dependent, title=title, init=init)
                          )

    def get_frame(self, index:int) -> FrameRecord:
        """
        Get frame.

        Args:
            index: Index of frame.
        Returns:
            Frame of record.
        Raises:
            ValueError:
                If index grate number of frames.
        """
        if index > len(self.frames) - 1:
            raise ValueError("Out of index")
        return self.frames[index]

    def clear(self) -> None:
        """
        Clear of the record.

        Returns:
            None
        """
        self.frames = []

    def get_size(self) -> int:
        """
        Get size of record.

        Returns:
            Number of frames.
        """
        return len(self.frames)

    def get_graph(self, index
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]|None:
        """
        Get graph by index.

        Args:
            index: Index of frame.
        Returns:
            Vertices, edges, and labels of edges.
        """
        while index > -1:
            frame = self.get_frame(index)
            graph = frame.get_graph()
            if not graph:
                index -= 1
                continue
            return graph
        return None

    def get_is_prev(self, index:int) -> bool:
        """
        Get the graph in the previous frame.

        Args:
            index: Index of frame.
        Returns:
            True if graph is None.
        """
        frame = self.get_frame(index)
        return frame.get_graph() is None

    def set_positions(self, positions:dict[str,np.array]) -> None:
        """
        Set positions.

        Args:
            positions: Positions of vertices.
        Returns:
            None
        """
        self.positions = positions

    def get_positions(self) -> dict[str,np.array]:
        """
        Get positions.

        Returns:
            Positions of vertices.
        """
        return self.positions

    def set_x_position_lighting(self, x_position_lighting:int) -> None:
        """
        Set x position of lighting.

        Args:
            x_position_lighting: x position of lighting.
        Returns:
            None
        """
        self.x_position_lighting = x_position_lighting

    def get_x_position_lighting(self) -> int:
        """
        Get x position of lighting.

        Returns:
            x position of lighting.
        """
        return self.x_position_lighting

def recording_graph(record:RecordGraph, collection:list[PauliString]=None,
                    lighting:PauliString=None, appending:PauliString=None,
                    contracting:PauliString=None, lits:list[PauliString]=None,
                    p:PauliString=None, q:PauliString = None,
                    removing_vertices:list[PauliString]=None,
                    replacing_vertices:list[PauliString]=None,
                    dependent:PauliString=None, title:str=None,
                    init:str=False) -> None:
    """
    Recording graph

    Args:
        record: Record of graph building.
        collection: List of vertices.
        lighting: Canonical graph join candidate.
        appending: Pauli string for appending.
        contracting: Pauli string for contracting.
        lits: List of lited vertices.
        p: Lited Pauli string in legs long 1.
        q: Unlited Pauli string in legs long 1.
        removing_vertices: List of removing vertices.
        replacing_vertices: List of replacing vertices.
        dependent: Dependent vertex.
        title: Title of frame.
        init: Init state of building.
    Returns:
        None
    """
    graph = None
    if collection is not None:
        vertices, edges, edge_labels = get_graph(collection)
        graph = FrameGraph(vertices, edges, edge_labels)
    record.append(graph, lighting=lighting, appending=appending,
                  contracting=contracting, lits=lits, p=p, q=q,
                  removing_vertices=removing_vertices,
                  replacing_vertices=replacing_vertices,
                  dependent=dependent, title=title, init=init)
