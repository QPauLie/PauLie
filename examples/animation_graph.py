from paulie.common.algebras import get_lie_algebra
from paulie.application.animation import animation_anti_commutation_graph


if __name__ == "__main__":
    animation_anti_commutation_graph(
        generators=get_lie_algebra()["a6"], 
        size=5
    )
