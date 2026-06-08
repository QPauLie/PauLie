import math
import os
import matplotlib.pyplot as plt
import networkx as nx
from paulie.helpers._recording import RecordGraph


def draw_single_frame(frame, ax, cmap, frame_idx):
    """Symulacja optyczna skanera: odchylenie refleksu czerwonego punktu."""
    G = frame.graph_state
    pos = nx.spring_layout(G, seed=42)

    if not pos:
        return

    # Srodek matrycy (punkt zero)
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    center_x = sum(xs) / len(xs)
    center_y = sum(ys) / len(ys)

    # Parametry eksperymentu: d = 1 metr, promien obrotu markera = 10cm
    dist_to_matrix = 1.0
    orbit_radius = 0.10

    # Ruch czerwonego przedmiotu po okregu w zaleznosci od klatki
    angle = (frame_idx * (2 * math.pi / 24))
    marker_x = center_x + orbit_radius * math.cos(angle)
    marker_y = center_y + orbit_radius * math.sin(angle)

    poly_nodes = []
    cross_nodes = []
    circle_nodes = []
    color_map_dict = {}

    # Obliczamy maksymalny odchyl, zeby dobrze wyskalowac barwy
    max_deviation = 0.0
    deviations = {}
    
    for node, p in pos.items():
        # Kat padania promienia z bableka na krecacy sie czerwony przedmiot
        dx = p[0] - marker_x
        dy = p[1] - marker_y
        # Odchylenie zalezy od odleglosci od osi prostopadlej i geometrii
        deviation = math.sqrt(dx**2 + dy**2) / dist_to_matrix
        deviations[node] = deviation
        if deviation > max_deviation:
            max_deviation = deviation

    # Mapowanie na automatyczna biblioteke kolorow (Zolty -> Czerwony -> Braz)
    for node in G.nodes():
        degree = G.degree(node)
        dev = deviations.get(node, 0.0)
        
        # Normowanie bledu (0.0 = zolty, 1.0 = ciemnoczerwony/braz)
        norm_dev = dev / max_deviation if max_deviation > 0 else 0.0
        color_map_dict[node] = cmap(norm_dev)

        if degree == 4:
            cross_nodes.append(node)
        elif degree > 4:
            poly_nodes.append(node)
        else:
            circle_nodes.append(node)

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#dddddd")

    # Rysujemy bable z automatycznym gradientem bledu optycznego
    if cross_nodes:
        nx.draw_networkx_nodes(
            G, pos, nodelist=cross_nodes, node_shape="d",
            node_size=750,
            node_color=[color_map_dict[n] for n in cross_nodes],
            ax=ax
        )
    if poly_nodes:
        nx.draw_networkx_nodes(
            G, pos, nodelist=poly_nodes, node_shape="p",
            node_size=700,
            node_color=[color_map_dict[n] for n in poly_nodes],
            ax=ax
        )
    if circle_nodes:
        nx.draw_networkx_nodes(
            G, pos, nodelist=circle_nodes, node_shape="o",
            node_size=600,
            node_color=[color_map_dict[n] for n in circle_nodes],
            ax=ax
        )

    nx.draw_networkx_labels(
        G, pos, ax=ax, font_size=10, font_color="black"
    )
    ax.set_title(f"Optyka skanera matrycy: Klatka {frame_idx:02d}")
    ax.axis("off")


def render_gradient_animation(record: RecordGraph, path: str):
    """Generuje sekwencje skanowania pod katem prostym z metra."""
    # Automatyczna biblioteka: Yellow-Orange-Red (Zolty -> Pomarancz -> Czerwony/Braz)
    cmap = plt.cm.YlOrRd
    
    base_dir = os.path.dirname(path)
    base_name = os.path.basename(path).replace(".gif", "").replace(".png", "")

    print("Uruchamiam skaner optyczny matrycy...")
    for idx, frame in enumerate(record.frames):
        fig, ax = plt.subplots(figsize=(6, 6))
        draw_single_frame(frame, ax, cmap, idx)
        
        frame_path = os.path.join(base_dir, f"{base_name}_{idx:02d}.png")
        fig.savefig(frame_path, dpi=100)
        plt.close()
        print(f"Wyznaczono odchylenie dla klatki: {idx:02d}")
