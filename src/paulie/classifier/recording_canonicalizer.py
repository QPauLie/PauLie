"""
    Recording wrapper around the canonicalizer.

    :class:`FrameRecorder` is the *concrete subscriber*: it reacts to events emitted by a
    :class:`~paulie.classifier.canonicalizer.Canonicalizer` and turns each of them into a frame in
    a :class:`~paulie.helpers._recording.RecordGraph`.

    :class:`RecordingCanonicalizer` is the *client*: it creates a canonicalizer and a recorder and
    wires them together. It does not change the algorithm, it only observes it, so it produces a
    :class:`~paulie.helpers._recording.RecordGraph` for any input the plain canonicalizer accepts.
"""
from __future__ import annotations

from typing import Any

from paulie.common.pauli_string_bitarray import PauliString
from paulie.classifier.canonicalizer import Canonicalizer
from paulie.classifier.observer import CanonicalizerObserver
from paulie.classifier.classification import Morph
from paulie.helpers._recording import RecordGraph, recording_graph


class FrameRecorder(CanonicalizerObserver):
    """
    Subscriber that records canonicalizer events as frames.
    """

    def __init__(self, record: RecordGraph) -> None:
        """
        Initialize the recorder.

        Args:
            record (RecordGraph): Record to append frames to.
        Returns:
            None
        """
        self.record = record

    def update(self, event: str, source: Canonicalizer, data: dict[str, Any]) -> None:
        """
        Translate an event into a frame.

        The current anticommutation graph is rebuilt from the canonicalizer state (or from an
        explicit ``collection`` for the initial frame). When a ``lighting`` vertex is present, the
        set of lit vertices is derived from the canonicalizer, so the recorder never needs the
        algorithm to compute them.

        Args:
            event (str): Name of the event, used as the frame title.
            source (Canonicalizer): The canonicalizer emitting the event.
            data (dict[str, Any]): Contextual data describing the event.
        Returns:
            None
        """
        collection = data.get("collection")
        if collection is None:
            collection = source._current_collection()  # pylint: disable=protected-access

        lighting = data.get("lighting")
        lits = None
        if lighting is not None:
            lits = [w for w in collection
                if w != lighting and source._is_lit(lighting, w)]  # pylint: disable=protected-access

        recording_graph(
            self.record,
            collection=collection,
            lighting=lighting,
            appending=data.get("appending"),
            contracting=data.get("contracting"),
            lits=lits,
            p=data.get("p"),
            q=data.get("q"),
            removing_vertices=self._as_list(data.get("removing")),
            replacing_vertices=self._as_list(data.get("replacing")),
            dependent=data.get("dependent"),
            title=event,
            init=data.get("init", False),
        )

    @staticmethod
    def _as_list(value: Any) -> list[PauliString] | None:
        """
        Normalize a role payload to a list of Pauli strings.

        Args:
            value (Any): A single Pauli string, a list of them, or None.
        Returns:
            list[PauliString] | None: List form of the payload, or None.
        """
        if value is None:
            return None
        if isinstance(value, list):
            return value
        return [value]


class RecordingCanonicalizer(Canonicalizer):
    """
    Canonicalizer that records its own transformation steps.

    Drop-in replacement for :class:`~paulie.classifier.canonicalizer.Canonicalizer`: it accepts the
    same input and returns the same :class:`~paulie.classifier.classification.Morph`, while emitting
    a frame for every step into a :class:`~paulie.helpers._recording.RecordGraph`.
    """

    def __init__(self, record: RecordGraph | None = None) -> None:
        """
        Initialize the recording canonicalizer.

        Args:
            record (RecordGraph, optional): Record to append frames to. Defaults to a fresh record.
        Returns:
            None
        """
        super().__init__()
        self.record = record if record is not None else RecordGraph()
        self.recorder = FrameRecorder(self.record)
        self.events.subscribe(self.recorder)

    def get_record(self) -> RecordGraph:
        """
        Get the recording produced by this canonicalizer.

        Returns:
            RecordGraph: Recording of the transformation steps.
        """
        return self.record

    def build_canonical_graph(self, vertex_stack: list[PauliString]) -> Morph:
        """
        Build a canonical graph from a stack of connected generators while recording.

        Args:
            vertex_stack (list[PauliString]): Generator stack.
        Returns:
            Morph: Canonical graph.
        """
        return super().build_canonical_graph(vertex_stack)
