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
        Constructor
        Args:
            vertices: list of vertices
            edges: list of edges
            edge_labels: list of edge's labels
        Returns: None
        """
        self.vertices = vertices
        self.edges = edges
        self.edge_labels = edge_labels

    def get_graph(self
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get graph pf the frame
        Args: empty
        Returns:
            the vertices, edges, and labels of edges
        """
        return self.vertices.copy(), self.edges.copy(), self.edge_labels

class FrameRecord:
    """
    Frame of record
    """
    def __init__(self, graph:FrameGraph=None, lighting:PauliString=None,
                 appending:PauliString=None,  contracting=PauliString,
                 lits:list[PauliString]=None, p:PauliString=None,
                 q:PauliString=None, removing_vertices:PauliString=None,
                 replacing_vertices:list[PauliString]=None, dependent:PauliString=None,
                 title:str=None, init:bool=False) -> None:
        """
        Constructor
        Args:
             graph: frame of graph
             lighting: canonical graph join candidate
             appending: Pauli string for appending
             contracting: Pauli string for contracting
             lits: list of lited vertices
             p: lited Pauli string in legs long 1
             q: unlited Pauli string in legs long 1
             removing_vertices: list of removing vertices
             replacing_vertices: list of replacing vertices
             dependent: dependent vertix
             title: title of frame
             init: init state of building
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
        Get graph of record's frame
        Args: empty
        Returns:
            the vertices, edges, and labels of edges
        """
        if not self.graph:
            return None
        return self.graph.get_graph()
    def get_lighting(self) -> str:
        """
        Get lighting
        Args: empty
        Returns:
           lighting
        """
        return self.lighting

    def get_title(self) -> str:
        """
        Get title
        Args: empty
        Returns:
            title
        """
        return self.title

    def is_appending(self) -> bool:
        """
        Check appending
        Args: empty
        Returns:
            True if not appending
        """
        return not self.appending

    def get_is_appending(self, vertix: str) -> bool:
        """
        Check vertix is appending
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is appending
        """
        if not self.appending:
            return False
        return self.appending == vertix

    def get_is_contracting(self, vertix:str) -> bool:
        """
        Check vertix is contracting
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is contracting

        """
        if not self.contracting:
            return False
        return self.contracting == vertix

    def get_is_p(self, vertix:str) -> bool:
        """
        Check vertix is p
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is lited in leg long 1
        """
        if not self.p:
            return False
        return self.p == vertix

    def get_is_q(self, vertix:str) -> bool:
        """
        Check vertix is q
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is unlited in leg long 1
        """
        if not self.q:
            return False
        return self.q == vertix

    def get_is_dependent(self, vertix:str) -> bool:
        """
        Check vertix is dependent
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is dependent
        """
        if not self.dependent:
            return False
        return self.dependent == vertix

    def get_is_lits(self, vertix:str) -> bool:
        """
        Check vertix is lited
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is lited
        """
        return vertix in self.lits

    def is_removing(self) -> bool:
        """
        Check is removing
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is not removing
        """
        return not self.removing_vertices

    def get_is_removing(self, vertix:str) -> bool:
        """
        Check vertix is removing
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is removing
        """
        return vertix in self.removing_vertices

    def get_is_replacing(self, vertix:str) -> bool:
        """
        Check vertix is replacing
        Args:
            vertix: vertix for checking
        Returns:
            True if vertix is replacing
        """
        return vertix in self.replacing_vertices

    def get_init(self) -> bool:
        """
        Get init frame 
        Args: empty
        Returns:
            True if init state
        """
        return self.init
class RecordGraph:
    """
    Record of graph
    """
    def __init__(self) -> None:
        """
        Constructor
        Args: empty
        Returns: None
        """
        self.frames: FrameRecord = []
        self.positions = {}
        self.x_position_lighting = 0


    def append_frame(self, frame: FrameRecord) -> None:
        """
        Append frame
        Args:
            frame: frame of record
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
        Make frame and append to record
        Args:
             graph: frame of graph
             lighting: canonical graph join candidate
             appending: Pauli string for appending
             contracting: Pauli string for contracting
             lits: list of lited vertices
             p: lited Pauli string in legs long 1
             q: unlited Pauli string in legs long 1
             removing_vertices: list of removing vertices
             replacing_vertices: list of replacing vertices
             dependent: dependent vertix
             title: title of frame
             init: init state of building
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
        Get frame
        Args:
            index: index of frame
        Returns:
            frame of record
        Raises:
              ValueError:
                 if index grate number of frames
        """
        if index > len(self.frames) - 1:
            raise ValueError("Out of index")
        return self.frames[index]

    def clear(self) -> None:
        """
        Clear of the record
        Args: empty
        Returns: None
        """
        self.frames = []

    def get_size(self) -> int:
        """
        Get size of record
        Args: empty
        Returns:
            number of frames
        """
        return len(self.frames)

    def get_graph(self, index
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]|None:
        """
        Get graph by index
        Args:
            index: index of frame
        Returns:
            the vertices, edges, and labels of edges
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
        Get the graph in the previous frame
        Args: 
            index: index of frame
        Returns:
            True if graph is None
        """
        frame = self.get_frame(index)
        return frame.get_graph() is None

    def set_positions(self, positions:dict[str,np.array]) -> None:
        """
        Set positions
        Args:
            positions: positions of vertices
        Returns:
            None
        """
        self.positions = positions

    def get_positions(self) -> dict[str,np.array]:
        """
        Get positions
        Args: empty
        Returns:
            positions of vertices
        """
        return self.positions

    def set_x_position_lighting(self, x_position_lighting:int) -> None:
        """
        Set x position of lighting
        Args:
            x_position_lighting: x position of lighting
        Returns:
            None
        """
        self.x_position_lighting = x_position_lighting

    def get_x_position_lighting(self) -> int:
        """
        Get x position of lighting
        Args: empty
        Returns:
            x position of lighting
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
        record: record of graph building
        collection: list of vertices
        lighting: canonical graph join candidate
        appending: Pauli string for appending
        contracting: Pauli string for contracting
        lits: list of lited vertices
        p: lited Pauli string in legs long 1
        q: unlited Pauli string in legs long 1
        removing_vertices: list of removing vertices
        replacing_vertices: list of replacing vertices
        dependent: dependent vertix
        title: title of frame
        init: init state of building
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
