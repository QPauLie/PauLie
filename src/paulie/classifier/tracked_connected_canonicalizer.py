from paulie.common.pauli_string_bitarray import PauliString
from paulie.classifier.connected_canonicalizer import ConnectedCanonicalizer
from paulie.classifier.classification import Morph

class TrackedConnectedCanonicalizer(ConnectedCanonicalizer):
    def __init__(self):
        self.multiplied_by: dict[PauliString, PauliString] = {}
        super().__init__()

    def _representative(self, v: PauliString):
        u = self.multiplied_by[v]
        if u is None:
            return v
        return u @ v

    def _tracked_multiply(self, v: PauliString, w: PauliString):
        u = self.multiplied_by[v]
        if u is None:
            u = w
        else:
            u = u @ w
        self.multiplied_by[v] = u
        return v

    def _dependency_check(self, length_1_legs: list[list[PauliString]]):
        # Hacky way to couple the representative information for undoing
        modified_legs = [[self._representative(leg[0]), self.multiplied_by[leg[0]]] for leg in length_1_legs]
        reduced_legs = super().dependency_check(modified_legs)
        independent_legs = []
        for leg in reduced_legs:
            if leg[1] is None:
                independent_legs.append([leg[0]])
            else:
                independent_legs.append([leg[1] @ leg[0]])
        return independent_legs

    def build_canonical_graph(self, vertex_stack: list[PauliString]) -> None:
        self.__init__()
        self._set_vertex_stack(vertex_stack)
        for v in vertex_stack:
            self.multiplied_by[v] = None
        super()._connected_canonical_graph(vertex_stack)

    def get_morph(self) -> Morph:
        """
        Get the canonical graph built by the canonicalizer.

        Returns:
            Morph: Canonical graph.
        """
        legs = self.legs.copy()
        legs.insert(0, [self.central_vertex])

        return Morph(legs, sum(legs, []), self.vertex_stack)
