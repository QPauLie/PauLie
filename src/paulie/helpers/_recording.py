"""
    Module that contains classes for creating a recording of the anticommutation graph
    transformation process.
"""
import numpy as np
from paulie.common.get_graph import get_graph
from paulie.common.pauli_string_bitarray import PauliString

class FrameGraph:
    """
    Stores the graph data in a frame.
    """
    def __init__(self, vertices:list[str], edges:list[tuple[str,str]],
                 edge_labels:dict[tuple[str,str],str] = None
    ) -> None:
        """
        Initialize the frame with graph data (vertices, edges, edge labels).

        Args:
            vertices (list[str]): List of vertices.
            edges (list[tuple[str,str]]): List of edges.
            edge_labels (dict[tuple[str,str],str]): List of edge labels.
        Returns:
            None
        """
        self.vertices = vertices
        self.edges = edges
        self.edge_labels = edge_labels

    def get_graph(self
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get the graph data in the frame.

        Returns:
            Vertices, edges, and labels of edges.
        """
        return self.vertices.copy(), self.edges.copy(), self.edge_labels

class FrameRecord:
    """
    A single frame of a recording.
    """
    def __init__(self, graph:FrameGraph=None, lighting:PauliString=None,
                 appending:PauliString=None,  contracting=PauliString,
                 lits:list[PauliString]=None, p:PauliString=None,
                 q:PauliString=None, removing_vertices:list[PauliString]=None,
                 replacing_vertices:list[PauliString]=None, dependent:PauliString=None,
                 title:str=None, init:bool=False) -> None:
        """
        Initialize the frame with graph data and the current state of the algorithm.

        Args:
            graph (FrameGraph, optional): Graph data.
            lighting (PauliString, optional): Vertex to be added which is generating the lighting.
            appending (PauliString, optional): Vertex to which the new vertex will be attached.
            contracting (PauliString): Vertex which is being contracted.
            lits (list[PauliString], optional): List of lit vertices.
            p (PauliString, optional): Lit vertex in a leg of length 1.
            q (PauliString, optional): Unlit vertex in a leg of length 1.
            removing_vertices (list[PauliString], optional): List of vertices to be removed.
            replacing_vertices (list[PauliString], optional): List of vertices to be replaced.
            dependent (PauliString, optional): Dependent vertex.
            title (str, optional): Title of frame.
            init (bool, optional): Whether this is the initial frame. Defaults to `False`.
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
        Get graph data stored in the frame.

        Returns:
            Vertices, edges, and labels of edges.
        """
        if not self.graph:
            return None
        return self.graph.get_graph()
    def get_lighting(self) -> str:
        """
        Get the vertex to be added.

        Returns:
            str: The vertex to be added.
        """
        return self.lighting

    def get_title(self) -> str:
        """
        Get the title of the frame.

        Returns:
            str: The title of the frame.
        """
        return self.title

    def is_appending(self) -> bool:
        """
        Check if the new vertex is being appended to the graph.

        Returns:
            bool: True if not appending.
        """
        return not self.appending

    def get_is_appending(self, vertex: str) -> bool:
        """
        Check if the vertex is being appended to.

        Args:
            vertex (str): Vertex for checking.
        Returns:
            bool: True if vertex is being appended to.
        """
        if not self.appending:
            return False
        return self.appending == vertex

    def get_is_contracting(self, vertex:str) -> bool:
        """
        Check if the vertex is being contracted.

        Args:
            vertex (str): Vertex for checking.
        Returns:
            bool: True if vertex is being contracted.

        """
        if not self.contracting:
            return False
        return self.contracting == vertex

    def get_is_p(self, vertex:str) -> bool:
        """
        Check if vertex is p.

        Args:
            vertex (str): Vertex for checking.
        Returns:
            bool: True if vertex is p.
        """
        if not self.p:
            return False
        return self.p == vertex

    def get_is_q(self, vertex:str) -> bool:
        """
        Check if vertex is q.

        Args:
            vertex (str): Vertex for checking.
        Returns:
            bool: True if vertex is q.
        """
        if not self.q:
            return False
        return self.q == vertex

    def get_is_dependent(self, vertex:str) -> bool:
        """
        Check if vertex is dependent.

        Args:
            vertex (str): Vertex for checking.
        Returns:
            bool: True if vertex is dependent.
        """
        if not self.dependent:
            return False
        return self.dependent == vertex

    def get_is_lits(self, vertex:str) -> bool:
        """
        Check if vertex is lit.

        Args:
            vertex (str): Vertex for checking.
        Returns:
            bool: True if vertex is lit.
        """
        return vertex in self.lits

    def is_removing(self) -> bool:
        """
        Check if vertices are not being removed.

        Returns:
            bool: True if vertices are not being removed.
        """
        return not self.removing_vertices

    def get_is_removing(self, vertex:str) -> bool:
        """
        Check if vertex is being removed.

        Args:
            vertex (str): Vertex for checking.
        Returns:
            bool: True if vertex is being removed.
        """
        return vertex in self.removing_vertices

    def get_is_replacing(self, vertex:str) -> bool:
        """
        Check if vertex is being replaced.

        Args:
            vertex (str): Vertex for checking.
        Returns:
            bool: True if vertex is being replaced.
        """
        return vertex in self.replacing_vertices

    def get_init(self) -> bool:
        """
        Check if this is the initial frame.

        Returns:
            bool: True if this is the initial frame.
        """
        return self.init

class RecordGraph:
    """
    Recording of the canonical graph transformation process.
    """
    def __init__(self) -> None:
        """
        Initialize an empty record.

        Returns:
            None
        """
        self.frames: FrameRecord = []
        self.positions = {}
        self.x_position_lighting = 0


    def append_frame(self, frame: FrameRecord) -> None:
        """
        Append a frame to the record.

        Args:
            frame (FrameRecord): Frame to be added to the record.
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
        Make a frame and append it to the record.
        Args:
            graph (FrameGraph, optional): Graph data.
            lighting (PauliString, optional): Vertex to be added which is generating the lighting.
            appending (PauliString, optional): Vertex to which the new vertex will be attached.
            contracting (PauliString, optional): Vertex which is being contracted.
            lits (list[PauliString], optional): List of lit vertices.
            p (PauliString, optional): Lit vertex in a leg of length 1.
            q (PauliString, optional): Unlit vertex in a leg of length 1.
            removing_vertices (list[PauliString], optional): List of vertices to be removed.
            replacing_vertices (list[PauliString], optional): List of vertices to be replaced.
            dependent (PauliString, optional): Dependent vertex.
            title (str, optional): Title of frame.
            init (bool, optional): Whether this is the initial frame. Defaults to `False`.
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
        Get a frame by its index.

        Args:
            index (int): Index of frame.
        Returns:
            FrameRecord: Frame of record.
        Raises:
            ValueError: If is index is greater than the number of frames.
        """
        if index > len(self.frames) - 1:
            raise ValueError("Out of index")
        return self.frames[index]

    def clear(self) -> None:
        """
        Clear the record.

        Returns:
            None
        """
        self.frames = []

    def get_size(self) -> int:
        """
        Get the number of frames in the record.

        Returns:
            int: Number of frames.
        """
        return len(self.frames)

    def get_graph(self, index
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]|None:
        """
        Get graph data by frame index.

        Args:
            index (int): Index of frame.
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
        Check if a frame has no graph data.

        Args:
            index (int): Index of frame.
        Returns:
            bool: True if the frame has no graph data.
        """
        frame = self.get_frame(index)
        return frame.get_graph() is None

    def set_positions(self, positions:dict[str,np.array]) -> None:
        """
        Set the positions of the vertices.

        Args:
            positions (dict[str,numpy.array]): Positions of vertices.
        Returns:
            None
        """
        self.positions = positions

    def get_positions(self) -> dict[str,np.array]:
        """
        Get the positions of the vertices.

        Returns:
            dict[str,numpy.array]: Positions of vertices.
        """
        return self.positions

    def set_x_position_lighting(self, x_position_lighting:int) -> None:
        """
        Set x position of the vertex to be added.

        Args:
            x_position_lighting (int): x position of the vertex to be added.
        Returns:
            None
        """
        self.x_position_lighting = x_position_lighting

    def get_x_position_lighting(self) -> int:
        """
        Get x position of the vertex to be added.

        Returns:
            x position of the vertex to be added.
        """
        return self.x_position_lighting

def recording_graph(record:RecordGraph, collection:list[PauliString]=None,
                    lighting:PauliString=None, appending:PauliString=None,
                    contracting:PauliString=None, lits:list[PauliString]=None,
                    p:PauliString=None, q:PauliString = None,
                    removing_vertices:list[PauliString]=None,
                    replacing_vertices:list[PauliString]=None,
                    dependent:PauliString=None, title:str=None,
                    init:bool=False) -> None:
    """
    Append a frame with the given data to the record.

    Args:
        record (RecordGraph): Record of graph building.
        collection (list[PauliString], optional): List of vertices.
        lighting (PauliString, optional): Vertex to be added which is generating the lighting.
        appending (PauliString, optional): Vertex to which the new vertex will be attached.
        contracting (PauliString, optional): Vertex which is being contracted.
        lits (list[PauliString], optional): List of lit vertices.
        p (PauliString, optional): Lit vertex in a leg of length 1.
        q (PauliString, optional): Unlit vertex in a leg of length 1.
        removing_vertices (list[PauliString], optional): List of vertices to be removed.
        replacing_vertices (list[PauliString], optional): List of vertices to be replaced.
        dependent (PauliString, optional): Dependent vertex.
        title (str, optional): Title of frame.
        init (bool, optional): Whether this is the initial frame. Defaults to `False`.
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
