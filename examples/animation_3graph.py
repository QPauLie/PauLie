from paulie.application.animation import animation_anti_commutation_graph


if __name__ == "__main__":
    animation_anti_commutation_graph(["XYI", "IXY", "XZY"], initGraph=True, storage={"filename":"data/example_a.html", "writer":"html"})
    animation_anti_commutation_graph(["XY", "XZ"], size = 4, initGraph=True, storage={"filename": "data/example_b.html", "writer": "html"})
