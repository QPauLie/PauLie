"""
Bulding animation
"""
from paulie import animation_anti_commutation_graph, get_pauli_string as p


if __name__ == "__main__":
    animation_anti_commutation_graph(p(["XYI", "IXY", "XZY"]),
                                     storage={"filename":"data/example_a.html", "writer":"html"})
    animation_anti_commutation_graph(p(["XY", "XZ"], n = 4),
                                     storage={"filename": "data/example_b.html", "writer": "html"})
