"""
Generate preview media for the classification animations.

Produces, for each worked example, an interactive HTML player, an animated GIF, and a still
PNG of every frame. Run with `uv run python docs/generate_animations.py`. The interactive
players embedded in the documentation are generated automatically at build time by
`docs/source/conf.py`; this script is a convenience for previewing the media locally.
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

from paulie import animation_anti_commutation_graph, get_pauli_string as p  # pylint: disable=wrong-import-position
from paulie.application.animation import _build_record  # pylint: disable=wrong-import-position

OUT = Path(__file__).resolve().parent / "preview"

EXAMPLES = {
    "example_c": p(["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]),
    "example_d": p(["XY", "XZ"], n=4),
}


def _slug(title: str) -> str:
    """
    Turn a frame title into a filename-safe slug.

    Args:
        title (str): Frame title.
    Returns:
        str: Filename-safe slug.
    """
    return re.sub(r"[^0-9a-zA-Z]+", "_", title).strip("_").lower()


def main() -> None:
    """
    Generate an HTML player, a GIF, and per-frame PNG stills for each example.

    Returns:
        None
    """
    OUT.mkdir(exist_ok=True)
    for name, generators in EXAMPLES.items():
        anim = animation_anti_commutation_graph(generators, interval=1200)
        # Interactive HTML player with play / pause / step / frame slider controls.
        (OUT / f"{name}.html").write_text(anim.to_jshtml(default_mode="loop"), encoding="utf-8")
        # Animated GIF.
        anim.save(str(OUT / f"{name}.gif"), writer="pillow", fps=1)
        # Still PNG of every frame, named by index and title.
        record = _build_record(generators)
        for i in range(record.get_size()):
            anim._func(i)  # pylint: disable=protected-access
            anim._fig.savefig(  # pylint: disable=protected-access
                OUT / f"{name}_{i:02d}_{_slug(record.get_frame(i).get_title())}.png", dpi=150)
        print(f"{name}: {record.get_size()} frames -> {OUT}")


if __name__ == "__main__":
    main()
