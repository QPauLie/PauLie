"""
Printing graph
"""
from paulie.common.pauli_string_bitarray import PauliString

def print_vertix(debug:bool, vertix:PauliString,
                 title:str="") -> None:
    """
    Printing vertix if debug.

    Args:
        debug: Debug flag.
        vertix: Vertix to print.
        title: Title to print.
    Returns:
        None
    """
    if debug:
        print(f"{title} {vertix}")

def print_vertices(debug:bool, vertices:list[PauliString],
                   title:str = "") -> None:
    """
    Printing list of vertices if debug.

    Args:
        debug: Debug flag.
        vertices: Vertices to print.
        title: Title to print.
    Returns:
        None
    """
    if debug is False:
        return

    print(f"----{title}--lenght = {len(vertices)}")
    for v in vertices:
        print_vertix(debug, v)
    print("-------------------")

def print_lit_vertices(debug:bool, vertices:list[PauliString],
                       lits:list[PauliString],
                       title:str = "") -> None:
    """
    Printing list of vertices with lits if debug.

    Args:
        debug: Debug flag.
        vertices: List of vertices to print.
        lits: List of lited vertices.
        title: Title to print.
    Returns:
        None
    """
    if debug is False:
        return
    print(f"----{title}--lenght = {len(vertices)}")
    for v in vertices:
        title = ""
        if v in lits:
            title = "*"
        print_vertix(debug, v, title)
    print("-------------------")


class Debug:
    """
    Debug class
    """
    def __init__(self, debug:bool) -> None:
        """
        Constructor.

        Args:
            debug: Debug flag.
        Returns:
            None
        """
        self.debug = debug
        self.save_debug = debug

    def get_debug(self) -> bool:
        """
        Get debug flag.

        Returns:
            True if debug flag is set in True.
        """
        return self.debug

    def set_debug(self, debug:bool) -> None:
        """
        Set debug flag.

        Args:
            debug: Debug flag.
        Returns:
            None
        """
        self.debug = debug

    def debuging(self) -> None:
        """
        Switch to debug mode.

        Returns: None
        """
        self.debug = True

    def restore(self) -> None:
        """
        Restore debug mode.

        Returns: None
        """
        self.debug = self.save_debug

    def print_vertix(self, vertix:PauliString,
                     title:str="") -> None:
        """
        Printing vertix if debug.

        Args:
            vertix: Vertix to print.
            title: Title to print.
        Returns:
            None
        """
        print_vertix(self.debug, vertix, title)


    def print_vertices(self, vertices: list[PauliString],
                       title:str="") -> None:
        """
        Printing list of vertices if debug.

        Args:
            vertices: List of vertices to print.
            title: Title to print.
        Returns:
            None
        """
        print_vertices(self.debug, vertices, title)

    def print_title(self, title:str) -> None:
        """
        Print title.

        Args:
            title: Title to print.
        Returns:
            None
        """
        if self.debug:
            if title != "":
                print(f"{title}")

    def print_lit_vertices(self,
                           vertices:list[PauliString],
                           lits:list[PauliString],
                           title:str = "") -> None:
        """
        Printing list of vertices with lits if debug.

        Args:
            vertices: List of vertices to print.
            lits: List of lited vertices.
            title: Title to print.
        Returns:
            None
        """
        print_lit_vertices(self.debug, vertices, lits, title)

    def is_pauli_string(self,
                        vertix:PauliString,
                        paulistring:PauliString) -> bool:
        """
        Pauli string equality check.

        Args:
            vertix: Vertix.
            paulistring: Pauli string to equality check.
        Returns:
           True if vertix equals Pauli string.
        """
        return paulistring == vertix
