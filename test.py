from paulie import get_pauli_string as p
from paulie.common.two_local_generators import G_LIE
from paulie.common.get_graph import get_graph
import matplotlib.pyplot as plt
import networkx as nx
import math

# Pick generator sets that yield canonical graphs of types A, B1, B2, B3
classes_to_test = {
    # "Type A": G_LIE["a22"],   # SO
    "Type B1": G_LIE["a9"],  # SP
    "Type B2": G_LIE["a11"], # SO (with long vertices)
    "Type B3": G_LIE["a12"]  # SU
}
# print(G_LIE["a1"])

for name, gen_list in classes_to_test.items():
    print(f"\n{'='*40}")
    print(f"Testing {name} with generators {gen_list}")
    print(f"{'='*40}")
    
    generators = p(gen_list, n=4)
    algebra = generators.get_algebra()
    if algebra is not None:
        print(f"Algebra: {algebra}")
    basis = generators.get_algebra_basis()
    if basis:
        print(f"Basis list length: {len(basis)}")
        print(f"Basis array shapes: {[b.shape for b in basis]}")
    else:
        print("Basis not implemented for this type.")
    
    vertices, edges, _ = generators.get_canonic_graph()
    print(f"Canonical Graph Vertices: {vertices}")
    print(f"Canonical Graph Edges: {edges}")

    # Plot the canonical graph with a star layout
    for morph in generators.get_class().morphs:
        one_legs, two_legs, _ = morph.counts()
        n_c = max(0, one_legs - 1)
        n_2 = two_legs
        print(f"Morph Stats -> n_c: {n_c}, n_2: {n_2}")
        
        if morph.is_empty():
            continue
        m_verts, m_edges, m_labels = get_graph(morph.get_vertices())
        graph = nx.Graph()
        graph.add_nodes_from(m_verts)
        graph.add_edges_from(m_edges)
        
        pos = {}
        if len(morph.legs) > 0:
            center_nodes = morph.legs[0]
            for j, node in enumerate(center_nodes):
                pos[str(node)] = (0.3 * j, 0)
                
            num_legs = len(morph.legs) - 1
            if num_legs > 0:
                angle_step = 2 * math.pi / num_legs
                for i, leg in enumerate(morph.legs[1:]):
                    angle = i * angle_step
                    for j, node in enumerate(leg):
                        r = j + 1
                        pos[str(node)] = (r * math.cos(angle), r * math.sin(angle))
                        
        for v in m_verts:
            if v not in pos:
                pos[v] = (0, 0)
                
        plt.figure(figsize=(6, 6))
        plt.title(f"{name}: {morph.get_algebra()}")
        nx.draw_networkx(
            graph, pos=pos, node_color='lightblue', 
            node_size=1500, font_size=9, font_weight='bold'
        )
        nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=m_labels, font_color="red")

plt.show()
