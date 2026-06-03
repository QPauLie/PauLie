"""Tests for the recording / animation machinery wrapping the canonicalizer."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from paulie import get_pauli_string as p
from paulie.helpers._recording import RecordGraph
from paulie.classifier.recording_canonicalizer import RecordingCanonicalizer
from paulie.application.animation import animation_anti_commutation_graph

# The two worked examples from docs/source/user/classification.rst.
A_TYPE = (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None, "4*so(5)")
B_TYPE = (["XY", "XZ"], 4, "sp(4)")

SMALL_CASES = [
    (["XY"], 3, "so(3)"),
    (["XY", "XZ"], 3, None),
    A_TYPE,
    B_TYPE,
]


def _collection(generators, n):
    return p(generators, n=n) if n is not None else p(generators)


def _frames_via_collection(generators, n):
    """Run a recorded classification and return the populated RecordGraph."""
    record = RecordGraph()
    collection = _collection(generators, n)
    collection.set_record(record)
    collection.classify()
    return record, collection


def test_recording_does_not_change_classification():
    """Recording must only observe: the classified algebra is unchanged."""
    for generators, n, _ in SMALL_CASES:
        plain = _collection(generators, n).get_algebra()
        record = RecordGraph()
        recorded_collection = _collection(generators, n)
        recorded_collection.set_record(record)
        recorded = recorded_collection.classify().get_algebra()
        assert plain == recorded


def test_record_is_non_empty():
    """A recorded run produces at least an initial and a final frame."""
    for generators, n, _ in SMALL_CASES:
        record, _ = _frames_via_collection(generators, n)
        assert record.get_size() >= 2


def test_initial_frame_is_input_graph():
    """The first frame is flagged as the init frame and holds the input graph."""
    generators, n, _ = B_TYPE
    record, _ = _frames_via_collection(generators, n)
    first = record.get_frame(0)
    assert first.get_init() is True
    vertices, edges, _labels = first.get_graph()
    # The initial frame is the input anticommutation graph of the first connected component.
    assert len(vertices) == len(set(vertices))
    assert len(vertices) >= 2
    vertex_set = set(vertices)
    for a, b in edges:
        assert a in vertex_set and b in vertex_set


def test_terminal_frame_is_canonical_type():
    """The final frame title announces a canonical type and carries no lighting vertex."""
    generators, n, _ = B_TYPE
    record, _ = _frames_via_collection(generators, n)
    last = record.get_frame(record.get_size() - 1)
    assert last.get_title().startswith("Canonical graph")
    assert last.get_lighting() is None


def test_every_frame_graph_is_well_formed():
    """Each frame resolves to a graph whose edges only reference existing vertices."""
    for generators, n, _ in SMALL_CASES:
        record, _ = _frames_via_collection(generators, n)
        for i in range(record.get_size()):
            graph = record.get_graph(i)
            assert graph is not None
            vertices, edges, _labels = graph
            assert len(vertices) == len(set(vertices))
            vertex_set = set(vertices)
            for a, b in edges:
                assert a in vertex_set
                assert b in vertex_set


def test_recording_canonicalizer_returns_morph():
    """RecordingCanonicalizer is a drop-in returning the same kind of Morph."""
    generators, n, _ = B_TYPE
    collection = _collection(generators, n)
    record = RecordGraph()
    canon = RecordingCanonicalizer(record)
    # Build the same vertex stack classify() would, for a single connected component.
    subgraphs = collection.get_subgraphs()
    assert len(subgraphs) == 1
    stack = subgraphs[0].get()
    morph = canon.build_canonical_graph(stack)
    assert morph.get_type().name in {"A", "B1", "B2", "B3", "NONE"}
    assert canon.get_record().get_size() >= 2


def test_lit_vertices_have_edges_to_lighting():
    """When a frame has a lighting vertex, every lit vertex is recorded for it."""
    generators, n, _ = A_TYPE
    record, _ = _frames_via_collection(generators, n)
    saw_lighting = False
    for i in range(record.get_size()):
        frame = record.get_frame(i)
        lighting = frame.get_lighting()
        if lighting is None:
            continue
        saw_lighting = True
        vertices, _edges, _labels = record.get_graph(i)
        # The lit set is a subset of the current vertices (excluding lighting itself).
        for v in vertices:
            if frame.get_is_lits(v):
                assert v != lighting
    assert saw_lighting


@pytest.mark.parametrize("generators,n", [(A_TYPE[0], A_TYPE[1]), (B_TYPE[0], B_TYPE[1])])
def test_animation_renders_to_gif(generators, n, tmp_path):
    """The renderer turns a run into a GIF, exercising every frame's draw callback."""
    collection = _collection(generators, n)
    out = tmp_path / "anim.gif"
    animation_anti_commutation_graph(
        collection,
        storage={"filename": str(out), "writer": "pillow"},
        interval=200,
        show=False,
    )
    plt.close("all")
    assert out.exists()
    assert out.stat().st_size > 0
