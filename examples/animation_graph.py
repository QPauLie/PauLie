from paulie.common.algebras import G_LIE
from paulie.application.animation import animation_anti_commutation_graph


if __name__ == "__main__":
    generators = G_LIE["a6"]
    animation_anti_commutation_graph(generators, size=5)
