"""
Printing graph
"""
from paulie.common.pauli_string_bitarray import PauliString

def print_vertix(debug:bool, vertix:PauliString,
                 title:str="") -> None:
    """
    Prnting vertix if debug
    """
    if debug:
        print(f"{title} {vertix}")

def print_vertices(debug:bool, vertices:list[PauliString],
                   title:str = "") -> None:
    """
    Prnting list of vertices if debug
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
    Prnting list of vertices with lits if debug
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
        Constuctor
        """
        self.debug = debug
        self.save_debug = debug

    def get_debug(self) -> bool:
        """
        Get debug flag
        """
        return self.debug

    def set_debug(self, debug:bool) -> None:
        """
        Set debug flag
        """
        self.debug = debug

    def debuging(self) -> None:
        """
        Switch to debug mode
        """
        self.debug = True

    def restore(self) -> None:
        """
        Restore debug mode
        """
        self.debug = self.save_debug

    def print_vertix(self, vertix:PauliString,
                     title:str="") -> None:
        """
        Prnting vertix if debug
        """
        print_vertix(self.debug, vertix, title)


    def print_vertices(self, vertices: list[PauliString],
                       title:str="") -> None:
        """
        Prnting list of vertices if debug
        """
        print_vertices(self.debug, vertices, title)

    def print_title(self, title:str) -> None:
        """
        Print title
        """
        if self.debug:
            if title != "":
                print(f"{title}")

    def print_lit_vertices(self,
                           vertices:list[PauliString],
                           lits:list[PauliString],
                           title:str = "") -> None:
        """
        Prnting list of vertices with lits if debug
        """
        print_lit_vertices(self.debug, vertices, lits, title)

    def is_pauli_string(self,
                        vertix:PauliString,
                        paulistring:PauliString) -> bool:
        """
        Pauli string equality check
        """
        return paulistring == vertix
