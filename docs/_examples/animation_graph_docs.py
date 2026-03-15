"""
Build animation files for the documentation.
"""

from pathlib import Path

from paulie import animation_anti_commutation_graph, get_pauli_string as p


ROOT = Path(__file__).resolve().parents[1]
MEDIA_DIR = ROOT / "docs" / "source" / "media"
MEDIA_DIR.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    anim = animation_anti_commutation_graph(
        p(["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]),
    )
    (MEDIA_DIR / "example_c.html").write_text(
        anim.to_jshtml(default_mode="once"),
        encoding="utf-8",
    )

    anim = animation_anti_commutation_graph(
        p(["XY", "XZ"], n=4),
    )
    (MEDIA_DIR / "example_d.html").write_text(
        anim.to_jshtml(default_mode="once"),
        encoding="utf-8",
    )
