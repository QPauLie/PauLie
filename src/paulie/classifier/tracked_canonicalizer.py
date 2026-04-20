"""
Tracked canonicalizer of generators
"""

from paulie.common.pauli_string_bitarray import PauliString
from paulie.classifier.canonicalizer import Canonicalizer
from paulie.classifier.classification import Morph

class TrackedCanonicalizer(Canonicalizer):
    """
    Class of tracked canonicalizer of generators
    """
    def __init__(self):
        """
        Initialize a TrackedCanonicalizer
        """
        self.multiplied_by: dict[PauliString, PauliString] = {}
        super().__init__()

    def _representative(self, v: PauliString) -> PauliString:
        """
        Pauli string representative
        Args:
            v (PauliString): Pauli string
        Returns:
            PauliString: representative
        """

        u = self.multiplied_by[v]
        if u is None:
            return v
        return u @ v

    def _tracked_multiply(self, v: PauliString, w: PauliString) -> PauliString:
        """
        Tracked multiply of two Pauli strings
        Args:
            v (PauliString): Pauli string
            w (PauliString): Pauli string
        Returns:
            PauliString: multiplication result
        """
        u = self.multiplied_by[v]
        if u is None:
            u = w
        else:
            u = u @ w
        self.multiplied_by[v] = u
        return v

    def _dependency_check(self, length_1_legs: list[list[PauliString]]):
        """
        Check depending on legs long 1
        Args:
            length_1_legs (list[list[PauliString]]): List of legs length 1
        Returns:
            list[list[PauliString]]: list of independent legs
        """
        # Hacky way to couple the representative information for undoing
        modified_legs = [[self._representative(leg[0]),
            self.multiplied_by[leg[0]]] for leg in length_1_legs]
        reduced_legs = super()._dependency_check(modified_legs)
        independent_legs = []
        for leg in reduced_legs:
            if leg[1] is None:
                independent_legs.append([leg[0]])
            else:
                independent_legs.append([leg[1] @ leg[0]])
        return independent_legs

    def _get_morph(self) -> Morph:
        """
        Get the canonical graph built by the canonicalizer.

        Returns:
            Morph: Canonical graph.
        """
        legs = self.legs.copy()
        legs.insert(0, [self.central_vertex])
        canonic_legs = [[self._representative(v) for v in leg] for leg in legs]
        return Morph(canonic_legs, sum(legs, []), self.vertex_stack)

    def build_canonical_graph(self, vertex_stack: list[PauliString]) -> Morph:
        """
        Build a canonical graph from a stack of connected generators

        Args:
            vertex_stack (list[PauliString]): Generator stack
        """
        self._set_vertex_stack(vertex_stack)
        for v in vertex_stack:
            self.multiplied_by[v] = None
        super()._connected_canonical_graph(vertex_stack)
        return self._get_morph()
