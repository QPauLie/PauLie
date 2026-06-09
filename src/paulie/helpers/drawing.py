"""
    Module with graph drawing utilities.
"""
import math
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.get_graph import get_graph
from paulie.helpers._recording import RecordGraph

#: Colour convention for the node roles tracked by the recorder.
#: The palette is a qualitative set chosen so the roles stay distinguishable from one another
#: (no two near-identical blues or magentas), which the docs legend renders as colour swatches.
NODE_ROLE_COLORS = {
    "lighting": "#e6194b",    # vertex currently being added (red)
    "dependent": "#9a6324",   # dependent vertex (brown)
    "contracting": "#f58231", # vertex currently being contracted (orange)
    "appending": "#3cb44b",   # attachment target / appended vertex (green)
    "removing": "#000000",    # vertex being temporarily removed (black)
    "replacing": "#911eb4",   # vertex being replaced (purple)
    "p": "#4363d8",           # lit vertex in a leg of length one (blue)
    "q": "#f032e6",           # unlit vertex in a leg of length one (pink)
    "lit": "#42d4f4",         # lit vertex (cyan)
    "unlit": "#d3d3d3",       # unlit vertex (light gray)
}

#: Human readable label for each node role, used to build the legend.
NODE_ROLE_LABELS = {
    "lighting": "vertex being added (lighting)",
    "dependent": "dependent vertex",
    "lit": "lit vertex",
    "unlit": "unlit vertex",
    "contracting": "vertex being contracted",
    "appending": "attachment target (appending)",
    "removing": "vertex being temporarily removed",
    "replacing": "vertex being replaced",
    "p": "lit vertex in a length-one leg (p)",
    "q": "unlit vertex in a length-one leg (q)",
}

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

def save_role_legend(filename: str) -> None:
    """
    Render the node-role colour legend as a standalone image.

    Each role is drawn as a filled colour swatch next to its description, so a reader can map the
    colours in an animation to their meaning. Colours come from :data:`NODE_ROLE_COLORS` and labels
    from :data:`NODE_ROLE_LABELS`.

    Args:
        filename (str): Path to write the legend image to.
    Returns:
        None
    """
    roles = list(NODE_ROLE_LABELS)
    fig, ax = plt.subplots(figsize=(5, 0.4 * len(roles) + 0.3))
    ax.axis("off")
    # One row per role, top to bottom: a square swatch on the left, its label to the right.
    for i, role in enumerate(roles):
        y = len(roles) - i
        ax.add_patch(plt.Rectangle((0, y - 0.4), 0.6, 0.6,
            facecolor=NODE_ROLE_COLORS[role], edgecolor="#666666"))
        ax.text(0.9, y - 0.1, NODE_ROLE_LABELS[role], va="center", fontsize=10)
    ax.set_xlim(0, 6)
    ax.set_ylim(0, len(roles) + 1)
    fig.tight_layout()
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)

def _node_color(frame, node: str, lighting: str | None) -> str:
    """
    Resolve the colour of a node for a frame from its role.

    Roles are checked in priority order so that the most specific role wins.

    Args:
        frame (FrameRecord): Frame being rendered.
        node (str): Vertex to colour.
        lighting (str | None): The lighting vertex of the frame, if any.
    Returns:
        str: A matplotlib colour.
    """
    if lighting is not None and node == lighting:
        if frame.get_is_dependent(node):
            return NODE_ROLE_COLORS["dependent"]
        return NODE_ROLE_COLORS["lighting"]
    if frame.get_is_removing(node):
        return NODE_ROLE_COLORS["removing"]
    if frame.get_is_dependent(node):
        return NODE_ROLE_COLORS["dependent"]
    if frame.get_is_replacing(node):
        return NODE_ROLE_COLORS["replacing"]
    if frame.get_is_appending(node):
        return NODE_ROLE_COLORS["appending"]
    if frame.get_is_contracting(node):
        return NODE_ROLE_COLORS["contracting"]
    if frame.get_is_p(node):
        return NODE_ROLE_COLORS["p"]
    if frame.get_is_q(node):
        return NODE_ROLE_COLORS["q"]
    if frame.get_is_lits(node):
        return NODE_ROLE_COLORS["lit"]
    return NODE_ROLE_COLORS["unlit"]

def _animation_graph(
    record: RecordGraph,
    interval: int = 1000,
    repeat: bool = False,
    storage: dict | None = None,
    show: bool = True,
) -> matplotlib.animation.Animation:
    """
    Animate the canonical graph construction from a recording.

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

    def clear() -> None:
        ax.clear()
        graph.remove_nodes_from(list(graph.nodes))

    def build_positions(
        edges: list[tuple[str, str]],
        center: str,
    ) -> tuple[dict[str, np.ndarray], float]:
        """
        Lay out a star-like canonical graph as 2D coordinates.

        The canonical graph is a central vertex with several legs (chains) hanging off it. The two
        longest legs are drawn as a single horizontal backbone through the centre, and the shorter
        legs fan out radially. This keeps the long leg readable even when it wraps onto several
        rows.

        Args:
            edges (list[tuple[str, str]]): Edges of the committed graph for this frame.
            center (str): The central vertex (legs are reconstructed by walking out from it).
        Returns:
            tuple[dict[str, np.ndarray], float]: Vertex coordinates and the x coordinate at which
            the incoming "lighting" vertex should be placed above the graph.
        """
        # Reconstruct the legs from the edge list. Each neighbour of the centre starts a leg...
        legs = []
        positions: dict[str, np.ndarray] = {}

        for edge in edges:
            if center in edge:
                v = edge[1] if center == edge[0] else edge[0]
                legs.append([v])

        # ...and we extend every leg outwards, hopping to the next unused neighbour until the
        # chain ends. This turns the edge set back into ordered chains of vertices.
        for leg in legs:
            current = leg[0]
            while True:
                is_found = False
                for edge in edges:
                    v = edge[1] if current == edge[0] else edge[0]
                    if current in edge and v not in leg and v != center:
                        leg.append(v)
                        current = v
                        is_found = True
                        break
                if not is_found:
                    break

        # Order legs shortest to longest so the two longest end up at the back of the list.
        legs.sort(key=len)
        n_legs_total = len(legs)

        max_line = 7          # vertices per row before the long leg wraps to a new row
        y_dist = 0.25         # vertical gap between wrapped rows
        pos_y = 0.0           # baseline y of the backbone
        y = pos_y
        center_x = 0.0
        x_position_lighting = 0.0
        x_first = 0.0
        x_last = 0.0

        # Horizontal spacing between adjacent vertices; tighter when there are many legs.
        # Must be defined for all branches.
        dist = 2.0 / (8 if n_legs_total > 7 else max(1, n_legs_total))

        if n_legs_total > 1:
            # Backbone: lay the two longest legs in one horizontal line through the centre.
            n = 0
            x = 1.0 + dist / 2.0

            # Second-longest leg fills the right half, placed right-to-left up to the centre.
            for v in reversed(legs[-2]):
                x -= dist
                if x_first == 0:
                    x_first = x
                positions[v] = np.array([x, y])
                n += 1

            # The centre sits between the two backbone legs.
            x -= dist
            positions[center] = np.array([x, y])
            center_x = x
            n += 1
            direction = -1

            # Longest leg continues to the left. When a row gets too long it wraps: drop down
            # one row and reverse horizontal direction so it snakes back (boustrophedon).
            for v in legs[-1]:
                if n > max_line:
                    if x_last == 0:
                        x_last = x
                    x += direction * dist
                    y -= y_dist
                    n = 0
                    max_line = 5
                    direction *= -1
                    if x_position_lighting == 0:
                        x_position_lighting = (x + 1 - dist / 2) / 2

                x += direction * dist
                positions[v] = np.array([x, y])
                n += 1
                if x_last == 0:
                    x_last = x

            # Park the incoming lighting vertex above the middle of the backbone.
            if x_position_lighting == 0:
                x_position_lighting = (x_first + x_last) / 2

            # The two longest legs are placed; the rest are handled radially below.
            legs = legs[:-2]

        if n_legs_total == 0:
            # Lone centre, nothing else to place.
            positions[center] = np.array([0.0, pos_y])
            x_position_lighting = 0.0

        elif n_legs_total == 1:
            # Single leg: just the centre and its one neighbour side by side.
            positions[legs[0][0]] = np.array([0.0, pos_y])
            positions[center] = np.array([0.25, pos_y])
            x_position_lighting = 0.125

        elif len(legs) > 0:
            # Remaining short legs fan out from the centre at evenly spaced angles, sweeping
            # back and forth so they spread above and below the backbone instead of overlapping.
            direction = 1
            n = len(legs)
            ang = 3 * math.pi / (2 * n)
            if ang >= math.pi / 2:
                ang = ang / 2
            c_ang = ang

            for leg in legs:
                # First vertex of the leg, one step out from the centre along the current angle.
                v = leg[0]
                y = pos_y + dist * math.sin(c_ang)
                x = center_x + dist * math.cos(c_ang)
                positions[v] = np.array([x, y])

                # A length-2 leg gets its second vertex one more step out along the same ray.
                if len(leg) > 1:
                    v = leg[1]
                    y = pos_y + 2 * dist * math.sin(c_ang)
                    x = center_x + 2 * dist * math.cos(c_ang)
                    positions[v] = np.array([x, y])

                # Advance the angle; once it passes the top, flip the sweep to the other side.
                c_ang += direction * ang
                if c_ang > 3 * math.pi / 4:
                    direction *= -1
                    c_ang = direction * ang

        return positions, x_position_lighting

    def update(num: int):
        clear()
        frame = record.get_frame(num)
        ax.set_title(f"{frame.get_title()}")

        vertices, edges, edge_labels = record.get_graph(num)
        vertices = list(vertices)
        edges = list(edges)

        center = None
        with_labels = True

        if len(vertices) > 0:
            center = vertices[0]
            if len(center) > 10:
                edge_labels = None
                with_labels = False

        lighting = frame.get_lighting()
        if lighting:
            if len(lighting) > 10:
                edge_labels = None
                with_labels = False

            if lighting not in vertices:
                vertices.append(lighting)

            for v in vertices:
                if frame.get_is_lits(v):
                    edges.append((lighting, v))

            if frame.get_is_dependent(lighting):
                dependent = lighting
                lighting = f"dependent {lighting}"
                edges.append((dependent, lighting))

        graph.add_nodes_from(vertices)
        graph.add_edges_from(edges)

        if frame.get_init() or center is None:
            positions = nx.spring_layout(graph)
        else:
            positions, x_position_lighting = build_positions(edges, center)
            if lighting:
                positions[lighting] = np.array([x_position_lighting, 1.0])

        # Robustness: make sure every drawn node has a position.
        missing = [node for node in graph if node not in positions]
        if missing:
            fallback = nx.spring_layout(graph, pos=positions, fixed=list(positions) or None)
            for node in missing:
                positions[node] = fallback[node]

        color_map = [_node_color(frame, node, lighting) for node in graph]

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

        return nx.draw_networkx(
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
