"""
    Regenerate the classification animations embedded in the documentation.

    Run from the repository root::

        python docs/generate_animations.py

    The animations are produced by the recording machinery
    (:func:`paulie.animation_anti_commutation_graph`) and written to ``docs/source/media``.
"""
import os

from paulie import get_pauli_string as p
from paulie import animation_anti_commutation_graph

MEDIA_DIR = os.path.join(os.path.dirname(__file__), "source", "media")

EXAMPLES = {
    # A-type canonical graph -> 4*so(5)
    "classification_a_type.gif": {
        "generators": ["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"],
        "n": None,
    },
    # B-type canonical graph (a_9) -> sp(4)
    "classification_b_type.gif": {
        "generators": ["XY", "XZ"],
        "n": 4,
    },
}


def main() -> None:
    """Generate every documentation animation."""
    for filename, spec in EXAMPLES.items():
        generators = (p(spec["generators"], n=spec["n"]) if spec["n"] is not None
                      else p(spec["generators"]))
        out = os.path.join(MEDIA_DIR, filename)
        print(f"Rendering {filename} (algebra = {generators.get_algebra()}) ...")
        animation_anti_commutation_graph(
            generators,
            storage={"filename": out, "writer": "pillow"},
            interval=1200,
            show=False,
        )
        print(f"  wrote {out}")


if __name__ == "__main__":
    main()
