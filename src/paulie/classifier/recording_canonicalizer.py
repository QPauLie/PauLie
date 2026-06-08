from typing import Any, Dict, List
from paulie.classifier.canonicalizer import Canonicalizer
from paulie.helpers._recording import RecordGraph


class RecordingCanonicalizer(Canonicalizer):
    """Szpieg zbierajacy snapshoty ewolucji bablekow."""

    def __init__(self):
        super().__init__()
        self.recorder = RecordGraph()
        self.current_roles: Dict[Any, str] = {}

    def record_step(self, event_name: str):
        """Zapisuje stan grafu tylko gdy on fizycznie istnieje."""
        if hasattr(self, "graph") and self.graph is not None:
            self.recorder.add_frame(
                event_name, self.graph, self.current_roles
            )

    def _set_vertex_stack(self, vertex_stack: list) -> None:
        """Inicjalizacja stosu - tu rodzi sie poczatkowy układ."""
        super()._set_vertex_stack(vertex_stack)
        # Jezeli baza stworzyla juz graf, robimy snapshot
        self.record_step("initial")

    def contract_vertex(self, vertex: Any):
        """Logika kurczenia bableka wzgledem srodka."""
        self.current_roles[vertex] = "contracting"
        self.record_step("before_contraction")

        super().contract_vertex(vertex)

        self.current_roles.pop(vertex, None)
        self.record_step("after_contraction")
