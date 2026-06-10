"""
Test recording and animation of the canonicalization process.
"""
import networkx as nx
import numpy as np
import pytest

from paulie import G_LIE, RecordingCanonicalizer, get_pauli_string as p
from paulie.application.animation import _build_record, animation_anti_commutation_graph
from paulie.classifier.classification import Classification
from paulie.helpers._recording import RecordGraph
from paulie.helpers.drawing import NODE_ROLE_COLORS, NODE_ROLE_LABELS, _stable_positions

A_TYPE = ["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]
B_TYPE = ["XY", "XZ"]

# (generators, qubit number, expected algebra)
SMALL_CASES = [
    (["XY"], 3, "so(3)"),
    (A_TYPE, None, "4*so(5)"),
    (B_TYPE, 4, "sp(4)"),
]


def _collection(generators: list[str], n: int | None):
    """
    Build a collection of generators for a test case.

    Args:
        generators (list[str]): Pauli strings.
        n (int): Number of qubits, or None.
    Returns:
        PauliStringCollection: Collection of generators.
    """
    return p(generators) if n is None else p(generators, n=n)


@pytest.mark.parametrize("generators, n, expected", SMALL_CASES)
def test_recording_matches_plain_classification(generators, n, expected) -> None:
    """
    Recording the construction does not change the classified algebra.
    """
    collection = _collection(generators, n)
    assert collection.get_algebra() == expected

    classification = Classification()
    record = RecordGraph()
    verts, edges, _ = collection.get_graph()
    graph = nx.Graph()
    graph.add_nodes_from(verts)
    graph.add_edges_from(edges)
    for component in nx.connected_components(graph):
        vertex_stack = [collection.create_instance(pauli_str=s)
            for s in nx.dfs_preorder_nodes(graph.subgraph(component))]
        vertex_stack.reverse()
        classification.add(RecordingCanonicalizer(record).build_canonical_graph(vertex_stack))
    assert str(classification.get_algebra()) == expected


@pytest.mark.parametrize("generators, n, _expected", SMALL_CASES)
def test_record_is_non_empty(generators, n, _expected) -> None:
    """
    A recording has at least an initial and a terminal frame.
    """
    record = _build_record(_collection(generators, n))
    assert record.get_size() >= 2


@pytest.mark.parametrize("generators, n, _expected", SMALL_CASES)
def test_initial_frame_is_input_graph(generators, n, _expected) -> None:
    """
    The initial frame is flagged as the input graph and is well-formed.
    """
    record = _build_record(_collection(generators, n))
    frame = record.get_frame(0)
    assert frame.get_init() is True
    vertices, edges, _ = frame.get_graph()
    assert len(set(vertices)) == len(vertices)
    assert len(vertices) >= 1
    for u, v in edges:
        assert u in vertices and v in vertices


@pytest.mark.parametrize("generators, n, expected", SMALL_CASES)
def test_terminal_frame_is_canonical(generators, n, expected) -> None:
    """
    The terminal frame is a canonical graph whose title carries the algebra name.
    """
    record = _build_record(_collection(generators, n))
    last = record.get_frame(record.get_size() - 1)
    assert last.get_title().startswith("Canonical graph")
    assert expected in last.get_title()
    assert last.get_lighting() is None


@pytest.mark.parametrize("generators, n, _expected", SMALL_CASES)
def test_every_frame_is_well_formed(generators, n, _expected) -> None:
    """
    Every frame resolves to a well-formed graph with unique vertices and valid edges.
    """
    record = _build_record(_collection(generators, n))
    for index in range(record.get_size()):
        vertices, edges, _ = record.get_graph(index)
        assert len(set(vertices)) == len(vertices)
        for u, v in edges:
            assert u in vertices and v in vertices


@pytest.mark.parametrize("generators, n, _expected", SMALL_CASES)
def test_lit_vertices_are_distinct_from_lighting(generators, n, _expected) -> None:
    """
    The vertex being added is never marked as one of its own lit neighbours.
    """
    record = _build_record(_collection(generators, n))
    for index in range(record.get_size()):
        frame = record.get_frame(index)
        lighting = frame.get_lighting()
        if lighting is not None:
            assert not frame.get_is_lits(lighting)


@pytest.mark.parametrize("generators, n, _expected", SMALL_CASES)
def test_recording_canonicalizer_is_drop_in(generators, n, _expected) -> None:
    """
    RecordingCanonicalizer returns a canonical Morph just like the base Canonicalizer.
    """
    collection = _collection(generators, n)
    verts, edges, _ = collection.get_graph()
    graph = nx.Graph()
    graph.add_nodes_from(verts)
    graph.add_edges_from(edges)
    for component in nx.connected_components(graph):
        vertex_stack = [collection.create_instance(pauli_str=s)
            for s in nx.dfs_preorder_nodes(graph.subgraph(component))]
        vertex_stack.reverse()
        morph = RecordingCanonicalizer().build_canonical_graph(vertex_stack)
        assert morph.get_type().name in {"A", "B1", "B2", "B3", "NONE"}


def test_layout_is_stable_when_graph_grows() -> None:
    """
    Placed vertices keep their position as the graph grows, so the animation does not jump.
    """
    positions: dict[str, np.ndarray] = {}
    angles: list[float] = []
    _stable_positions(positions, angles, [("C", "A"), ("C", "B")], "C")
    snapshot = {v: pos.copy() for v, pos in positions.items()}
    # Extend a leg and add a new leg.
    _stable_positions(positions, angles, [("C", "A"), ("A", "D"), ("C", "B"), ("C", "E")], "C")
    for v, pos in snapshot.items():
        assert np.allclose(positions[v], pos)
    assert "D" in positions and "E" in positions


def test_layout_is_stable_under_relabel() -> None:
    """
    A relabelled leg vertex inherits the position of the slot it occupies.
    """
    positions: dict[str, np.ndarray] = {}
    angles: list[float] = []
    _stable_positions(positions, angles, [("C", "A"), ("A", "B")], "C")
    pos_b = positions["B"].copy()
    # `B` is relabelled to `B2` at the same depth on the same leg.
    _stable_positions(positions, angles, [("C", "A"), ("A", "B2")], "C")
    assert np.allclose(positions["B2"], pos_b)


def test_layout_is_stable_across_a_fold() -> None:
    """
    Growing a leg past the fold length never moves the vertices already placed on it.
    """
    positions: dict[str, np.ndarray] = {}
    angles: list[float] = []
    chain = [f"V{i}" for i in range(14)]
    snapshots = {}
    for length in range(1, len(chain) + 1):
        edges = [("C", chain[0])] + [(chain[i], chain[i + 1]) for i in range(length - 1)]
        _stable_positions(positions, angles, edges, "C")
        for v, pos in snapshots.items():
            assert np.allclose(positions[v], pos)
        snapshots = {v: pos.copy() for v, pos in positions.items()}
    # The leg folds into rows instead of running off along a single ray.
    xs = [abs(pos[0]) for pos in positions.values()]
    assert max(xs) <= 8 * 0.6 + 1e-9


def test_role_legend_is_consistent() -> None:
    """
    Every node role has both a colour and a label, and the colours are hex values.
    """
    assert set(NODE_ROLE_LABELS) == set(NODE_ROLE_COLORS)
    for color in NODE_ROLE_COLORS.values():
        assert color.startswith("#") and len(color) == 7


@pytest.mark.slow
@pytest.mark.parametrize("name", sorted(G_LIE))
def test_recording_matches_plain_for_all_two_local_algebras(name) -> None:
    """
    Recording never changes the classified algebra, for every two-local algebra.
    """
    collection = p(G_LIE[name], n=6)
    expected = collection.get_algebra()
    classification = Classification()
    verts, edges, _ = collection.get_graph()
    graph = nx.Graph()
    graph.add_nodes_from(verts)
    graph.add_edges_from(edges)
    for component in nx.connected_components(graph):
        vertex_stack = [collection.create_instance(pauli_str=s)
            for s in nx.dfs_preorder_nodes(graph.subgraph(component))]
        vertex_stack.reverse()
        classification.add(RecordingCanonicalizer().build_canonical_graph(vertex_stack))
    assert sorted(str(classification.get_algebra()).split("+")) == sorted(expected.split("+"))


@pytest.mark.slow
@pytest.mark.parametrize("generators, n", [(A_TYPE, None), (B_TYPE, 4)])
def test_animation_renders(generators, n) -> None:
    """
    The animation renders to an interactive HTML player without error.
    """
    anim = animation_anti_commutation_graph(_collection(generators, n))
    html = anim.to_jshtml()
    assert isinstance(html, str) and len(html) > 0
