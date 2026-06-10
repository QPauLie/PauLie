"""
    Module with graph drawing utilities.
"""
import math
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
from paulie.helpers._recording import RecordGraph
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.get_graph import get_graph

# Colour of each node role in the animation. The same hex values are used in the colour
# legend of docs/source/user/classification.rst. The palette is chosen to be distinguishable
# (Sasha Trubetskoy's qualitative set) so the roles can be told apart at a glance.
NODE_ROLE_COLORS = {
    "lighting": "#e6194b",
    "dependent": "#9a6324",
    "contracting": "#f58231",
    "appending": "#3cb44b",
    "removing": "#000000",
    "replacing": "#911eb4",
    "p": "#4363d8",
    "q": "#f032e6",
    "lit": "#42d4f4",
    "unlit": "#d3d3d3",
}

# Human-readable description of each node role, used in the colour legend.
NODE_ROLE_LABELS = {
    "lighting": "Vertex being added (lighting)",
    "dependent": "Dependent vertex",
    "contracting": "Vertex being contracted",
    "appending": "Attachment target",
    "removing": "Vertex being removed",
    "replacing": "Vertex being replaced",
    "p": "Lit vertex in a leg of length one",
    "q": "Unlit vertex in a leg of length one",
    "lit": "Lit vertex",
    "unlit": "Unlit vertex",
}

# Constant spacing between consecutive vertices on a leg, and the fixed height at which the
# vertex being added hovers above the centre. Keeping these constant (rather than rescaling to
# fit the figure) is what stops vertices and edge lengths from jumping between frames.
_LAYOUT_DIST = 0.6
_LIGHTING_Y = 1.6
# Long legs snake into rows of at most _FOLD_LENGTH vertices, _FOLD_DIST apart. The fold is a
# pure function of the depth on the leg, so it is part of the layout from the first frame and
# never repositions a vertex when the leg grows.
_FOLD_LENGTH = 8
_FOLD_DIST = 0.5
# Angular slots for legs, in the order they are first assigned: the first two legs form a
# horizontal backbone through the centre, the rest fan out above and below.
_LEG_ANGLES = [
    0.0, math.pi,
    math.pi / 2, -math.pi / 2,
    math.pi / 4, -math.pi / 4, 3 * math.pi / 4, -3 * math.pi / 4,
    math.pi / 6, -math.pi / 6, 5 * math.pi / 6, -5 * math.pi / 6,
    math.pi / 3, -math.pi / 3, 2 * math.pi / 3, -2 * math.pi / 3,
]


def _reconstruct_legs(edges: list[tuple[str, str]], center: str) -> list[list[str]]:
    """
    Reconstruct the legs of a starlike graph from its edge list.

    Each neighbour of the centre starts a leg, which is extended by following the unique chain
    of edges away from the centre until it ends.

    Args:
        edges (list[tuple[str, str]]): Edges of the graph.
        center (str): Central vertex.
    Returns:
        list[list[str]]: Legs, each ordered from the centre outward.
    """
    legs = []
    for edge in edges:
        if center in edge:
            v = edge[1] if center == edge[0] else edge[0]
            legs.append([v])
    for leg in legs:
        current = leg[0]
        while True:
            is_found = False
            for edge in edges:
                if current in edge:
                    v = edge[1] if current == edge[0] else edge[0]
                    if v not in leg and v != center:
                        leg.append(v)
                        current = v
                        is_found = True
                        break
            if not is_found:
                break
    return legs


def _next_leg_angle(assigned_angles: list[float]) -> float:
    """
    Pick the next unused angular slot for a new leg.

    Args:
        assigned_angles (list[float]): Angles already assigned to legs.
    Returns:
        float: Angle for the new leg, in radians.
    """
    for angle in _LEG_ANGLES:
        if all(abs(angle - used) > 1e-6 for used in assigned_angles):
            return angle
    # Fall back to the golden angle once the fixed slots are exhausted.
    return (len(assigned_angles) * 2.399963) % (2 * math.pi)


def _leg_offset(depth: int, ray: np.ndarray) -> np.ndarray:
    """
    Offset of the vertex at a given depth on a leg, relative to the centre.

    The leg snakes: rows of ``_FOLD_LENGTH`` vertices run along the leg's ray, alternating
    direction, each row shifted perpendicular to the ray. The offset depends only on the depth,
    so a vertex keeps its position no matter when the leg grows past a fold.

    Args:
        depth (int): Index of the vertex on the leg, counted from the centre.
        ray (numpy.ndarray): Unit direction of the leg.
    Returns:
        numpy.ndarray: Offset from the centre position.
    """
    row, col = divmod(depth, _FOLD_LENGTH)
    along = col + 1 if row % 2 == 0 else _FOLD_LENGTH - col
    perp = np.array([-ray[1], ray[0]])
    return ray * along * _LAYOUT_DIST - perp * row * _FOLD_DIST


def _recover_leg_ray(leg: list[str], positions: dict[str, np.ndarray],
    center_pos: np.ndarray) -> np.ndarray | None:
    """
    Recover the outward unit direction of a leg from any already-placed vertex on it.

    Inverts :func:`_leg_offset` for the vertex's depth: the offset is linear in the ray, so the
    ray is the solution of a 2x2 system. A relabelled or newly added vertex on an existing leg
    is then positioned consistently from any surviving vertex, so placed vertices never move.

    Args:
        leg (list[str]): Leg vertices, ordered from the centre outward.
        positions (dict[str, numpy.ndarray]): Known vertex positions.
        center_pos (numpy.ndarray): Position of the central vertex.
    Returns:
        numpy.ndarray|None: Unit direction of the leg, or None if no vertex is placed yet.
    """
    for depth, v in enumerate(leg):
        if v in positions:
            row, col = divmod(depth, _FOLD_LENGTH)
            along = col + 1 if row % 2 == 0 else _FOLD_LENGTH - col
            # offset = along * d * ray - row * f * perp(ray) is linear in the ray:
            # offset = (along * d * I - row * f * J) @ ray with J the 90-degree rotation.
            a = along * _LAYOUT_DIST
            b = row * _FOLD_DIST
            matrix = np.array([[a, b], [-b, a]])
            ray = np.linalg.solve(matrix, positions[v] - center_pos)
            norm = float(np.linalg.norm(ray))
            if norm > 1e-9:
                return ray / norm
    return None


def _remember_angle(angle: float, assigned_angles: list[float]) -> None:
    """
    Record an angle as used, unless an equivalent angle is already recorded.

    Args:
        angle (float): Angle in radians.
        assigned_angles (list[float]): Angles already assigned to legs, mutated in place.
    Returns:
        None
    """
    for used in assigned_angles:
        if abs((angle - used + math.pi) % (2 * math.pi) - math.pi) < 1e-6:
            return
    assigned_angles.append(angle)


def _stable_positions(positions: dict[str, np.ndarray], assigned_angles: list[float],
    edges: list[tuple[str, str]], center: str) -> float:
    """
    Lay out the canonical graph, moving a vertex only when its place in the graph has changed.

    The centre sits at the origin, each leg lies along a fixed ray, and a vertex sits at the
    offset :func:`_leg_offset` of its depth on the leg, snaking into rows once the leg grows
    long. A leg keeps the ray recovered from its already-placed vertices, so growing a leg or
    relabelling a vertex in place never moves the others; only when the algorithm restructures
    a leg (splitting it or moving vertices between legs) are the affected vertices repositioned
    onto their new slots. Spacing is constant throughout, so distances never jump.

    Args:
        positions (dict[str, numpy.ndarray]): Persistent vertex positions, mutated in place.
        assigned_angles (list[float]): Angles already assigned to legs, mutated in place.
        edges (list[tuple[str, str]]): Edges of the current canonical graph.
        center (str): Central vertex.
    Returns:
        float: x coordinate of the central vertex.
    """
    if center not in positions:
        positions[center] = np.array([0.0, 0.0])
    center_pos = positions[center]

    claimed: list[np.ndarray] = []
    for leg in _reconstruct_legs(edges, center):
        ray = _recover_leg_ray(leg, positions, center_pos)
        # Two legs cannot share a ray: a split leg whose vertices still sit on the old ray
        # is given a fresh slot instead of overlapping the leg that kept it.
        if ray is not None and any(float(np.dot(ray, c)) > 0.999 for c in claimed):
            ray = None
        if ray is None:
            angle = _next_leg_angle(assigned_angles)
            ray = np.array([math.cos(angle), math.sin(angle)])
        _remember_angle(math.atan2(ray[1], ray[0]), assigned_angles)
        claimed.append(ray)
        for depth, v in enumerate(leg):
            expected = center_pos + _leg_offset(depth, ray)
            if v not in positions or not np.allclose(positions[v], expected, atol=1e-9):
                positions[v] = expected
    return float(center_pos[0])


def _view_limits(positions: dict[str, np.ndarray]
) -> tuple[tuple[float, float], tuple[float, float]] | None:
    """
    Compute fixed axis limits covering every vertex position, with a margin.

    Fixing the viewport across frames keeps the graph stationary so it grows in place instead
    of the view jumping or zooming between frames.

    Args:
        positions (dict[str, numpy.ndarray]): All vertex positions in the animation.
    Returns:
        tuple|None: ``((xmin, xmax), (ymin, ymax))`` limits, or None if there are no positions.
    """
    if not positions:
        return None
    xs = [pos[0] for pos in positions.values()]
    ys = [pos[1] for pos in positions.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys + [_LIGHTING_Y + _LAYOUT_DIST])
    margin_x = 0.15 * (xmax - xmin) + 0.2
    margin_y = 0.15 * (ymax - ymin) + 0.2
    return (xmin - margin_x, xmax + margin_x), (ymin - margin_y, ymax + margin_y)


def plot_graph(vertices:list[str],
               edges:list[tuple[str,str]],
               edge_labels:dict[tuple[str,str],str] = None) -> None:
    """
    Plot an arbitrary graph.

    Args:
        vertices (list[str]): List of vertices.
        edges (list[tuple[str, str]]): List of edges.
        edge_labels (dict[tuple[str,str],str], optional): List of edge labels. Defaults to None, in
            which case no edge labels are drawn.
    Returns:
        None
    """
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    graph.add_edges_from(edges)
    pos = nx.spring_layout(graph)
    if edge_labels is not None:
        nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=edge_labels,font_color="red")
    nx.draw_networkx(graph, pos=pos)
    #plt.show()

def plot_graph_by_nodes(nodes:PauliStringCollection,
                        commutators:PauliStringCollection|list[PauliString]=None
) -> None:
    """
    Plot anticommutator graph with edges labeled by commutator of endpoints.

    Args:
        nodes (PauliStringCollection): Generating set for anticommutator graph.
        commutators (PauliStringCollection|list[PauliString], optional): Only show those edges which
            have labels belonging to this set. Defaults to None, in which case all edges are shown.
    Returns:
        None
    """
    if not commutators:
        commutators = []
    vertices, edges, edge_labels = get_graph(nodes, commutators)
    return plot_graph(vertices, edges, edge_labels)

def _animation_graph(
    record: RecordGraph,
    interval: int = 1000,
    repeat: bool = False,
    storage: dict | None = None,
    show: bool = True,
) -> matplotlib.animation.Animation:
    """
    Animate the canonical graph construction from a recording.

    Vertex positions are cached in the recording and assigned incrementally, so a vertex keeps
    its position for the whole animation and the graph grows without jumping between frames.

    Args:
        record (RecordGraph): A recording of the canonical graph construction.
        interval (int, optional): Interval between frames in milliseconds.
        repeat (bool, optional): Whether to loop the animation.
        storage (dict, optional): Output location and format. Expected keys:
            - "filename": path to the output file
            - "writer": matplotlib writer name or writer instance
        show (bool, optional): Whether to display the animation window.

    Returns:
        matplotlib.animation.Animation
    """
    graph = nx.Graph()
    fig, ax = plt.subplots(figsize=(6, 4))

    # Pre-compute every frame's positions so the viewport can be fixed: the graph then grows
    # in place instead of the view jumping or zooming between frames. The layout state is
    # reset afterwards so playback rebuilds the identical position sequence from scratch.
    all_positions: dict[str, np.ndarray] = {}
    for num in range(record.get_size()):
        if record.get_frame(num).get_init():
            continue
        frame_vertices, frame_edges, _ = record.get_graph(num)
        if frame_vertices:
            _stable_positions(record.get_positions(), record.assigned_angles,
                              list(frame_edges), frame_vertices[0])
            all_positions.update(record.get_positions())
    view_limits = _view_limits(all_positions)
    record.set_positions({})
    record.assigned_angles.clear()

    def clear() -> None:
        ax.clear()
        graph.remove_nodes_from(list(graph.nodes))

    def update(num: int):
        clear()
        frame = record.get_frame(num)
        ax.set_title(f"{frame.get_title()}")

        vertices, base_edges, edge_labels = record.get_graph(num)
        vertices = list(vertices)
        base_edges = list(base_edges)

        center = vertices[0] if vertices else None
        # Labels become an unreadable smear on long Pauli strings or crowded graphs; the node
        # colours and the frame title still tell the story, so drop the labels there.
        with_labels = True
        if (center is not None and len(center) > 10) or len(vertices) > 12:
            edge_labels = None
            with_labels = False

        # The initial frame is the raw input graph, drawn with a spring layout. Every later
        # frame uses the persistent, incremental layout so placed vertices never move.
        if frame.get_init() or center is None:
            graph.add_nodes_from(vertices)
            graph.add_edges_from(base_edges)
            positions = nx.spring_layout(graph, seed=7)
            center_x = 0.0
        else:
            center_x = _stable_positions(record.get_positions(), record.assigned_angles,
                                         base_edges, center)
            positions = dict(record.get_positions())
            graph.add_nodes_from(vertices)
            graph.add_edges_from(base_edges)

        # Overlay the vertex being added at a fixed position above the centre, with edges to
        # the vertices it lights up.
        lighting = frame.get_lighting()
        if lighting:
            if len(lighting) > 10:
                edge_labels = None
                with_labels = False
            graph.add_node(lighting)
            positions[lighting] = np.array([center_x, _LIGHTING_Y])
            for v in vertices:
                if frame.get_is_lits(v):
                    graph.add_edge(lighting, v)
            if frame.get_is_dependent(lighting):
                dependent = lighting
                lighting = f"dependent {lighting}"
                graph.add_node(lighting)
                positions[lighting] = np.array([center_x, _LIGHTING_Y + _LAYOUT_DIST])
                graph.add_edge(dependent, lighting)

        if frame.is_appending() and frame.is_removing():
            graph.remove_nodes_from([v for v in vertices if frame.get_is_removing(v)])

        # Any node still without a position (a relabelled straggler) is relaxed in while every
        # already-placed vertex stays fixed.
        missing = [n for n in graph.nodes if n not in positions]
        if missing:
            fixed = [n for n in positions if n in graph]
            positions = nx.spring_layout(graph, pos=positions, fixed=fixed or None, seed=7)

        color_map = []
        for node in graph:
            if lighting == node:
                if frame.get_is_dependent(node):
                    color_map.append(NODE_ROLE_COLORS["dependent"])
                else:
                    color_map.append(NODE_ROLE_COLORS["lighting"])
            elif not frame.get_is_lits(node):
                if frame.get_is_q(node):
                    color_map.append(NODE_ROLE_COLORS["q"])
                elif frame.get_is_removing(node):
                    color_map.append(NODE_ROLE_COLORS["removing"])
                else:
                    color_map.append(NODE_ROLE_COLORS["unlit"])
            else:
                if frame.get_is_appending(node):
                    color_map.append(NODE_ROLE_COLORS["appending"])
                elif frame.get_is_contracting(node):
                    color_map.append(NODE_ROLE_COLORS["contracting"])
                elif frame.get_is_p(node):
                    color_map.append(NODE_ROLE_COLORS["p"])
                elif frame.get_is_replacing(node):
                    color_map.append(NODE_ROLE_COLORS["replacing"])
                else:
                    color_map.append(NODE_ROLE_COLORS["lit"])

        ax.axis("off")

        if edge_labels is not None:
            nx.draw_networkx_edge_labels(
                graph,
                pos=positions,
                edge_labels=edge_labels,
                font_color="red",
                hide_ticks=True,
                node_size=60,
                font_size=6,
                ax=ax,
            )

        artists = nx.draw_networkx(
            graph,
            pos=positions,
            node_color=color_map,
            hide_ticks=True,
            node_size=60,
            font_size=6,
            with_labels=with_labels,
            edge_color="#aaaaaa",
            ax=ax,
        )

        # Fix the viewport on the construction frames so the graph stays put as it grows.
        if not frame.get_init() and view_limits is not None:
            ax.set_xlim(view_limits[0])
            ax.set_ylim(view_limits[1])

        return artists

    ani = matplotlib.animation.FuncAnimation(
        fig,
        update,
        frames=record.get_size(),
        interval=interval,
        repeat=repeat,
    )

    if storage is not None:
        ani.save(filename=storage["filename"], writer=storage["writer"])

    if show:
        plt.show()

    return ani
