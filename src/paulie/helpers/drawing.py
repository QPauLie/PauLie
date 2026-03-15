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
        legs = []
        positions: dict[str, np.ndarray] = {}

        for edge in edges:
            if center in edge:
                v = edge[1] if center == edge[0] else edge[0]
                legs.append([v])

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

        legs.sort(key=len)
        n_legs_total = len(legs)

        max_line = 7
        y_dist = 0.25
        pos_y = 0.0
        y = pos_y
        center_x = 0.0
        x_position_lighting = 0.0
        x_first = 0.0
        x_last = 0.0

        # Must be defined for all branches.
        dist = 2.0 / (8 if n_legs_total > 7 else max(1, n_legs_total))

        if n_legs_total > 1:
            n = 0
            x = 1.0 + dist / 2.0

            for v in reversed(legs[-2]):
                x -= dist
                if x_first == 0:
                    x_first = x
                positions[v] = np.array([x, y])
                n += 1

            x -= dist
            positions[center] = np.array([x, y])
            center_x = x
            n += 1
            direction = -1

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

            if x_position_lighting == 0:
                x_position_lighting = (x_first + x_last) / 2

            # Remove the two longest legs; they are already placed.
            legs = legs[:-2]

        if n_legs_total == 0:
            positions[center] = np.array([0.0, pos_y])
            x_position_lighting = 0.0

        elif n_legs_total == 1:
            positions[legs[0][0]] = np.array([0.0, pos_y])
            positions[center] = np.array([0.25, pos_y])
            x_position_lighting = 0.125

        elif len(legs) > 0:
            direction = 1
            n = len(legs)
            ang = 3 * math.pi / (2 * n)
            if ang >= math.pi / 2:
                ang = ang / 2
            c_ang = ang

            for leg in legs:
                v = leg[0]
                y = pos_y + dist * math.sin(c_ang)
                x = center_x + dist * math.cos(c_ang)
                positions[v] = np.array([x, y])

                if len(leg) > 1:
                    v = leg[1]
                    y = pos_y + 2 * dist * math.sin(c_ang)
                    x = center_x + 2 * dist * math.cos(c_ang)
                    positions[v] = np.array([x, y])

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

            vertices.append(lighting)

            for v in vertices:
                if frame.get_is_lits(v):
                    edges.append((lighting, v))

            if frame.get_is_dependent(lighting):
                dependent = lighting
                lighting = f"dependent {lighting}"
                edges.append((dependent, lighting))

        if frame.is_appending() and frame.is_removing():
            vertices = [v for v in vertices if not frame.get_is_removing(v)]

        graph.add_nodes_from(vertices)
        graph.add_edges_from(edges)

        if frame.get_init():
            positions = nx.spring_layout(graph)
        else:
            if not record.get_is_prev(num):
                positions, x_position_lighting = build_positions(edges, center)
                record.set_positions(positions)
                record.set_x_position_lighting(x_position_lighting)
                if lighting:
                    positions[lighting] = np.array([x_position_lighting, 1.0])
            else:
                positions = record.get_positions()
                x_position_lighting = record.get_x_position_lighting()
                if lighting:
                    positions[lighting] = np.array([x_position_lighting, 1.0])

        color_map = []
        for node in graph:
            if lighting == node:
                if frame.get_is_dependent(node):
                    color_map.append("#2F4F4F")
                else:
                    color_map.append("red")
            else:
                if not frame.get_is_lits(node):
                    if frame.get_is_q(node):
                        color_map.append("#FF00FF")
                    else:
                        if frame.get_is_removing(node):
                            color_map.append("black")
                        else:
                            color_map.append("#cccccc")
                else:
                    if frame.get_is_appending(node):
                        color_map.append("#00FF00")
                    elif frame.get_is_contracting(node):
                        color_map.append("#008080")
                    else:
                        if frame.get_is_p(node):
                            color_map.append("#6A5ACD")
                        else:
                            if frame.get_is_replacing(node):
                                color_map.append("#8B008B")
                            else:
                                color_map.append("cyan")

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
