"""
    Module for animating the transformation of the anticommutation graph into a canonical form.
"""
import networkx as nx
from matplotlib.animation import Animation
from paulie.helpers._recording import RecordGraph
from paulie.helpers.drawing import _animation_graph
from paulie.classifier.recording_canonicalizer import RecordingCanonicalizer
from paulie.common.pauli_string_collection import PauliStringCollection


def _build_record(generators: PauliStringCollection) -> RecordGraph:
    """
    Record the canonicalization of a generating set.

    The anticommutation graph is split into its connected components and each component is
    canonicalized in turn, recording every transformation step into a single recording. This
    mirrors the component iteration of :func:`PauliStringCollection.classify`.

    Args:
        generators (PauliStringCollection): Collection of Pauli strings.
    Returns:
        RecordGraph: Recording of the canonical graph construction.
    """
    record = RecordGraph()
    verts, edges, _ = generators.get_graph()
    graph: nx.Graph[str] = nx.Graph()
    graph.add_nodes_from(verts)
    graph.add_edges_from(edges)
    for component in nx.connected_components(graph):
        vertex_stack = [generators.create_instance(pauli_str=s)
            for s in nx.dfs_preorder_nodes(graph.subgraph(component))]
        vertex_stack.reverse()
        RecordingCanonicalizer(record).build_canonical_graph(vertex_stack)
    return record


def animation_anti_commutation_graph(
    generators: PauliStringCollection,
    storage: dict[str, str] | None = None,
    interval: int = 1000,
    show: bool = False,
) -> Animation:
    """
    Animate the transformation of the anticommutation graph into a canonical form.

    Args:
        generators (PauliStringCollection): Collection of Pauli strings.
        storage (dict[str, str], optional): Output location and format. Expected keys
            ``filename`` and ``writer``. Defaults to None, in which case the animation is not
            saved.
        interval (int, optional): Interval between frames in milliseconds.
        show (bool, optional): Whether to display the animation window.
    Returns:
        matplotlib.animation.Animation
    """
    record = _build_record(generators)
    return _animation_graph(
        record,
        interval=interval,
        storage=storage,
        show=show,
    )
