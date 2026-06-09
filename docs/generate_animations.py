"""
    Regenerate the classification animations embedded in the documentation.

    Run from the repository root::

        python docs/generate_animations.py

    For each example this writes three assets into ``docs/source/media``:

    - a ``.gif`` for a quick non-interactive preview,
    - a ``.html`` interactive player (matplotlib ``to_jshtml``) with play / pause / step / loop
      controls, which is what the docs embed,

    and once, a ``classification_legend.png`` colour legend shared by both examples.
"""
import os

from paulie import get_pauli_string as p
from paulie import animation_anti_commutation_graph
from paulie.helpers.drawing import save_role_legend

MEDIA_DIR = os.path.join(os.path.dirname(__file__), "source", "media")

EXAMPLES = {
    # A-type canonical graph -> 4*so(5)
    "classification_a_type": {
        "generators": ["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"],
        "n": None,
    },
    # B-type canonical graph (a_9) -> sp(4)
    "classification_b_type": {
        "generators": ["XY", "XZ"],
        "n": 4,
    },
}


def main() -> None:
    """Generate the legend and every documentation animation."""
    legend_path = os.path.join(MEDIA_DIR, "classification_legend.png")
    save_role_legend(legend_path)
    print(f"wrote {legend_path}")

    for name, spec in EXAMPLES.items():
        generators = (p(spec["generators"], n=spec["n"]) if spec["n"] is not None
                      else p(spec["generators"]))
        print(f"Rendering {name} (algebra = {generators.get_algebra()}) ...")

        # Build the animation once, then export it as both a gif and an interactive player.
        ani = animation_anti_commutation_graph(generators, interval=1200, show=False)

        gif_path = os.path.join(MEDIA_DIR, f"{name}.gif")
        ani.save(filename=gif_path, writer="pillow")
        print(f"  wrote {gif_path}")

        html_path = os.path.join(MEDIA_DIR, f"{name}.html")
        with open(html_path, "w", encoding="utf-8") as handle:
            handle.write(ani.to_jshtml(default_mode="once"))
        print(f"  wrote {html_path}")


if __name__ == "__main__":
    main()
