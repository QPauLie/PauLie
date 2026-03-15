"""
Building animation.
"""

from paulie import animation_anti_commutation_graph, get_pauli_string as p


if __name__ == "__main__":
    animation_anti_commutation_graph(
        p(["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]),
        storage={
            "filename": "./docs/source/media/example_c.html",
            "writer": "html",
        },
    )
    animation_anti_commutation_graph(
        p(["XY", "XZ"], n=4),
        storage={
            "filename": "./docs/source/media/example_d.html",
            "writer": "html",
        },
    )