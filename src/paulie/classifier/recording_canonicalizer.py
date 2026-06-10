"""
Recording canonicalizer of generators.
"""
from paulie.common.pauli_string_bitarray import PauliString
from paulie.classifier.canonicalizer import Canonicalizer
from paulie.classifier.classification import Classification, Morph
from paulie.classifier.observer import CanonicalizerEvent, CanonicalizerObserver
from paulie.helpers._recording import RecordGraph, recording_graph


class _FrameRecorder(CanonicalizerObserver):
    """
    Observer that records each canonicalization event as a frame.
    """
    def __init__(self, record: RecordGraph) -> None:
        """
        Initialize a frame recorder.

        Args:
            record (RecordGraph): Recording to append frames to.
        Returns:
            None
        """
        self.record = record

    def update(self, event: CanonicalizerEvent, data: dict) -> None:
        """
        Record a canonicalization event as a frame.

        Args:
            event (CanonicalizerEvent): Type of event emitted.
            data (dict): Contextual data associated with the event.
        Returns:
            None
        """
        title = data.get("title")
        if title is None:
            lighting = data.get("lighting")
            title = f"{event.value}: {lighting}" if lighting is not None else event.value
        recording_graph(self.record, collection=data.get("collection"),
                        lighting=data.get("lighting"), lits=data.get("lits"),
                        p=data.get("p"), q=data.get("q"),
                        removing_vertices=data.get("removing"),
                        dependent=data.get("dependent"), title=title,
                        init=bool(data.get("init", False)))


class RecordingCanonicalizer(Canonicalizer):
    """
    Canonicalizer that records the transformation process into a RecordGraph.
    """
    def __init__(self, record: RecordGraph | None = None) -> None:
        """
        Initialize a RecordingCanonicalizer.

        Args:
            record (RecordGraph, optional): Recording to append frames to. Defaults to None, in
                which case a new recording is created.
        Returns:
            None
        """
        super().__init__()
        self.record = record if record is not None else RecordGraph()
        self.events.subscribe(_FrameRecorder(self.record))

    def get_record(self) -> RecordGraph:
        """
        Get the recording of the transformation process.

        Returns:
            RecordGraph: Recording of the canonical graph construction.
        """
        return self.record

    def build_canonical_graph(self, vertex_stack: list[PauliString]) -> Morph:
        """
        Build a canonical graph from a stack of connected generators, recording each step.

        Args:
            vertex_stack (list[PauliString]): Generator stack.
        Returns:
            Morph: Canonical graph.
        """
        self._notify(CanonicalizerEvent.ANTICOMMUTATION_GRAPH,
                     collection=list(vertex_stack), init=True)
        morph = super().build_canonical_graph(vertex_stack)
        self._notify(CanonicalizerEvent.CANONICAL,
                     collection=self._current_vertices(), title=self._final_title(morph))
        return morph

    def _final_title(self, morph: Morph) -> str:
        """
        Build the title of the terminal frame including the classified algebra.

        Args:
            morph (Morph): Canonical graph.
        Returns:
            str: Terminal frame title of the form `Canonical graph of type B3: su(1024)`.
        """
        classification = Classification()
        classification.add(morph)
        return (f"Canonical graph of type {morph.get_type().name}: "
                f"{classification.get_algebra()}")
